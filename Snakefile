# Snakefile for mapping reads from Whole Genome NGS

# All input parameters are being read from config.yaml
# Script will work on files in "project_dir"
# genome in FASTA format is located in genome_dir and is prefixed with
# "genome_prefix"like:
# It assumes reads are located in "reads_dir" and they have the structure:
# {sample}/{file_name}_{pair}.fq.gz, where {pair} may be "1" or "2"

project_main = config["project_dir"]
project_reads = "/".join((project_main, config["reads_dir"]))
project_genome = "/".join((project_main, config["genome_dir"], config["genome_prefix"]))
project_variants = "/".join((project_main, config["variants_dir"], config["genome_prefix"]))
project_samples = "/".join((project_main, config["samples_dir"]))
reads_pairs = ["1", "2"]

genome_prefix=config["genome_prefix"]
variants_dir = "/".join((project_main, config["variants_dir"]))
snpeff_database = config["snpeff_database"]

project_logs = "/".join((project_main, "logs"))

(picard_version) = glob_wildcards("/opt/conda/share/picard-{version}/picard.jar")
if not len(picard_version.version):
    raise ValueError("no Picard jar file found")
picard_jar = "/opt/conda/share/picard-{}/picard.jar".format(picard_version.version[0])
print("Using Picard version {}".format(picard_jar))

(snpeff_version) = glob_wildcards("/opt/conda/share/snpeff-{version}/snpEff.jar")
if not len(snpeff_version.version):
    raise ValueError("no SnpEff jar file found")
snpeff_jar = "/opt/conda/share/snpeff-{}/snpEff.jar".format(picard_version.version[0])
print("Using SnpEff version {}".format(snpeff_version)

def read_samples_names():
    (dirs, ) = glob_wildcards(project_samples + "/{dir}")
    dirs = [ele for ele in dirs if "/" not in ele]
    return dirs

samples_names = read_samples_names()

print("{} samples loaded".format(len(samples_names)))
print(*samples_names, sep=", ")

threads_max = config["threads_max"]
memory_max = config["memory_max"]

rule fastqc:
    input:
        ancient("{project_samples}/{sample}/{sample}_{pair}.fq.gz")
    output:
        "{project_samples}/{sample}/metrics/fastqc/{sample}_{pair}/{sample}_{pair}_fastqc.zip",
    threads: 1
    params:
        output_dir="{project_samples}/{sample}/metrics/fastqc/{sample}_{pair}"
    shell:
        "fastqc {input} -o {params.output_dir}"

rule prepare_genome_index:
    input:
        ancient(project_genome + ".fa")
    output:
        protected(project_genome + ".fa.fai")
    run:
        "samtools faidx {input} -o {output}"

rule prepare_genome_dictionary:
    input:
        fa=ancient(project_genome + ".fa"),
        fai=rules.prepare_genome_index.output
    output:
        protected(project_genome + ".dict")
    params:
        picard_jar = picard_jar
    threads: 1
    run:
        "java -jar {params.picard_jar} CreateSequenceDictionary \
        R={input} O={output.dict}"

rule prepare_bwa_genome:
    input:
        ancient(project_genome + ".fa")
    output:
        amb=protected(project_genome + ".amb"),
        ann=protected(project_genome + ".ann"),
        bwt=protected(project_genome + ".bwt"),
        pac=protected(project_genome + ".pac"),
        sa=protected(project_genome + ".sa")
    threads: 1
    shell:
        "bwa index -p " + project_genome + " {input}"

rule bwa_map_reads:
    priority: 2000
    input:
        reads_1=ancient("{project_samples}/{sample}/{sample}_1.fq.gz"),
        reads_2=ancient("{project_samples}/{sample}/{sample}_2.fq.gz"),
        genome=ancient(rules.prepare_bwa_genome.output)
    output:
        temp("{project_samples}/{sample}/{sample}.unsorted.bam")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.bwa.log")
    params:
        prefix=project_genome,
        rg=r"@RG\tID:{sample}\tSM:{sample}\tPL:illumina\tLB:{sample}\tPU:{sample}"
    threads: threads_max
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {params.prefix} \
        {input.reads_1} {input.reads_2} | samtools view -Sb - > {output}) 2> {log}"

rule samtools_sort_bwa_map:
    priority: 1950
    input:
        rules.bwa_map_reads.output
    output:
        temp("{project_samples}/{sample}/{sample}.sorted.bam")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.sort_bwa_map.log")
    threads: 10
    resources:
        mem_mb=4096
    shell:
        "samtools sort {input} \
        -o {output} \
        -@ {threads} \
        -m {resources.mem_mb}M \
        -T /tmp/{wildcards.sample} 2> {log}"

rule samtools_index_sorted_bwa_map:
    priority: 1900
    input:
        rules.samtools_sort_bwa_map.output
    output:
        temp("{project_samples}/{sample}/{sample}.sorted.bam.bai")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.index_sorted_bwa_map.log")
    threads: 4
    shell:
        "samtools index -@ {threads} {input} {output}"

rule picard_mark_duplicates:
    priority: 1850
    input:
        sorted_bam=rules.samtools_sort_bwa_map.output,
        sorted_bai=rules.samtools_index_sorted_bwa_map.output,
    output:
        marked_bam=protected("{project_samples}/{sample}/{sample}.bam"),
        marked_metrics=protected("{project_samples}/{sample}/metrics/{sample}.MarkDuplicates.txt")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.MarkDuplicates.log")
    params:
        picard_jar=picard_jar
    threads: threads_max
    resources:
        mem_mb = memory_max
    shell:
        "java -Xmx{resources.mem_mb}m -jar {params.picard_jar} \
         MarkDuplicates \
         I={input.sorted_bam} \
         O={output.marked_bam} \
         M={output.marked_metrics} \
         2> {log}"

rule samtools_index_marked:
    priority: 1800
    input:
        rules.picard_mark_duplicates.output.marked_bam
    output:
        protected("{project_samples}/{sample}/{sample}.bam.bai")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.index_marked.log")
    threads: 4
    shell:
        "samtools index -@ {threads} {input} {output}"

rule samtools_flagstats:
    priority: 1500
    input:
        marked_bam=rules.picard_mark_duplicates.output.marked_bam,
        marked_bai=rules.samtools_index_marked.output
    output:
        protected("{project_samples}/{sample}/metrics/{sample}.flagstats")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.flagstats.log")
    threads: threads_max
    shell:
        "(samtools flagstat -@ {threads} {input.marked_bam} > {output}) 2> {log}"

rule samtools_idxstats:
    priority: 1500
    input:
        marked_bam=rules.picard_mark_duplicates.output.marked_bam,
        marked_bai=rules.samtools_index_marked.output
    output:
        protected("{project_samples}/{sample}/metrics/{sample}.idxstats")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.idxstats.log")
    threads: 1
    shell:
        "(samtools idxstats {input.marked_bam} > {output}) 2> {log}"

rule samtools_stats:
    priority: 1500
    input:
        marked_bam=rules.picard_mark_duplicates.output.marked_bam,
        marked_bai=rules.samtools_index_marked.output
    output:
        protected("{project_samples}/{sample}/metrics/{sample}.stats")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.stats.log")
    threads: 1
    shell:
        "(samtools stats {input.marked_bam} > {output}) 2> {log}"

rule picard_size_metrics:
    priority: 1500
    input:
        marked_bam=rules.picard_mark_duplicates.output.marked_bam,
        marked_bai=rules.samtools_index_marked.output
    output:
        txt=protected("{project_samples}/{sample}/metrics/{sample}.InsertSizeMetrics.txt"),
        pdf=protected("{project_samples}/{sample}/metrics/{sample}.InsertSizeMetrics.pdf")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.InsertSizeMetrics.log")
    params:
        picard_jar = picard_jar
    threads: 1
    resources:
        mem_mb = 8192
    shell:
        "java -Xmx{resources.mem_mb}m -jar {params.picard_jar} \
        CollectInsertSizeMetrics \
        I={input.marked_bam} \
        O={output.txt} \
        H={output.pdf} \
        M=0.5 2> {log}"

rule picard_alignment_summary:
    priority: 1500
    input:
        marked_bam=rules.picard_mark_duplicates.output.marked_bam,
        marked_bai=rules.samtools_index_marked.output,
        genome=ancient(project_genome + ".fa"),
        genome_dict=rules.prepare_genome_dictionary.output
    output:
        txt=protected("{project_samples}/{sample}/metrics/{sample}.AlignmentSummaryMetrics.txt")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.AlignmentSummaryMetrics.log")
    params:
        picard_jar = picard_jar
    threads: 1
    resources:
        mem_mb = 8192
    shell:
        "java -Xmx{resources.mem_mb}m -jar {params.picard_jar} \
        CollectAlignmentSummaryMetrics \
        R={input.genome} \
        I={input.marked_bam} \
        O={output.txt} 2> {log}"

rule picard_wgs_metrics:
    priority: 1500
    input:
        marked_bam=rules.picard_mark_duplicates.output.marked_bam,
        marked_bai=rules.samtools_index_marked.output,
        genome=ancient(project_genome + ".fa"),
        genome_dict=rules.prepare_genome_dictionary.output
    output:
        txt=protected("{project_samples}/{sample}/metrics/{sample}.WgsMetrics.txt")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.WgsMetrics.log")
    params:
        picard_jar = picard_jar
    threads: 4
    resources:
        mem_mb = 32768
    shell:
        "java -Xmx{resources.mem_mb}m -jar {params.picard_jar} \
        CollectWgsMetrics \
        R={input.genome} \
        I={input.marked_bam} \
        O={output.txt} 2> {log}"

bait_bed="/".join((project_main, config["bait_intervals"] + ".bed"))
bait_intervals="/".join((project_main, config["bait_intervals"] + ".intervals"))

rule prepare_bait_intervals:
    input: 
        bed=bait_bed,
        genome_dict=rules.prepare_genome_dictionary.output
    output: bait_intervals
    log: project_logs + "/Picard.BaitBedToIntervalList.log"
    params:
        picard_jar=picard_jar
    threads: 1
    resources: mem_mb=8192
    shell:
        "java -Xmx{resources.mem_mb}m -jar {params.picard_jar} \
        BedToIntervalList \
        I={input.bed} \
        SD={input.genome_dict} \
        O={output} \
        2> {log}"

target_bed="/".join((project_main, config["target_intervals"] + ".bed"))
target_intervals="/".join((project_main, config["target_intervals"] + ".intervals"))

rule prepare_target_intervals:
    input: 
        bed=target_bed,
        genome_dict=rules.prepare_genome_dictionary.output
    output: target_intervals
    log: project_logs + "/Picard.TargetBedToIntervalList.log"
    params:
        picard_jar=picard_jar
    threads: 1
    resources: mem_mb=8192
    shell:
        "java -Xmx{resources.mem_mb}m -jar {params.picard_jar} \
        BedToIntervalList \
        I={input.bed} \
        SD={input.genome_dict} \
        O={output} \
        2> {log}"

rule picard_hs_metrics:
    priority: 1500
    input:
        marked_bam=rules.picard_mark_duplicates.output.marked_bam,
        marked_bai=rules.samtools_index_marked.output,
        genome=ancient(project_genome + ".fa"),
        genome_dict=rules.prepare_genome_dictionary.output,
        target_intervals=rules.prepare_target_intervals.output,
        bait_intervals=rules.prepare_bait_intervals.output
    output:
        txt=protected("{project_samples}/{sample}/metrics/{sample}.HsMetrics.txt")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.HsMetrics.log")
    params:
         picard_jar=picard_jar
    threads:2 
    resources:
         mem_mb=16384
    shell:
          "java -Xmx{resources.mem_mb}m -jar {params.picard_jar} \
          CollectHsMetrics \
          R={input.genome} \
          I={input.marked_bam} \
          TARGET_INTERVALS={input.target_intervals} \
          BAIT_INTERVALS={input.bait_intervals} \
          O={output.txt}\
          2> {log}"

rule picard_gc_bias_metrics:
    priority: 1500
    input:
        marked_bam=rules.picard_mark_duplicates.output.marked_bam,
        marked_bai=rules.samtools_index_marked.output,
        genome=ancient(project_genome + ".fa"),
        genome_dict=rules.prepare_genome_dictionary.output
    output:
        metrics=protected("{project_samples}/{sample}/metrics/{sample}.GCBiasMetrics.txt"),
        chart=protected("{project_samples}/{sample}/metrics/{sample}.GCBiasMetrics.pdf"),
        summary=protected("{project_samples}/{sample}/metrics/{sample}.GCBiasMetricsSummary.txt")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.GCBiasMetrics.log")
    params:
        picard_jar = picard_jar
    threads: 1
    resources:
        mem_mb = 8192
    shell:
        "java -Xmx{resources.mem_mb}m -jar {params.picard_jar} \
        CollectGcBiasMetrics \
        R={input.genome} I={input.marked_bam} O={output.metrics} S={output.summary} \
        CHART={output.chart} 2> {log}"

rule picard_validate_sam_file:
    priority: 1500
    input:
        marked_bam=rules.picard_mark_duplicates.output.marked_bam,
        marked_bai=rules.samtools_index_marked.output,
        genome=ancient(project_genome + ".fa"),
        genome_dict=rules.prepare_genome_dictionary.output
    output:
        txt=protected("{project_samples}/{sample}/metrics/{sample}.ValidateSamFile.txt")
    log:
        protected("{project_samples}/{sample}/metrics/{sample}.ValidateSamFile.log")
    params:
        picard_jar = picard_jar
    threads: 1
    resources:
        mem_mb = 8192
    shell:
        "java -Xmx{resources.mem_mb}m -jar {params.picard_jar} \
        ValidateSamFile \
        I={input.marked_bam} OUTPUT={output.txt} \
        R={input.genome} MODE='VERBOSE' 2> {log}"

rule gatk_index_variants:
    input:
        "{}.vcf.gz".format(project_variants)
    output:
        "{}.vcf.gz.tbi".format(project_variants)
    log:
        "{}.log".format(project_variants)
    shell:
        "gatk IndexFeatureFile \
        -F {input} \
        -O {output} 2> {log}"

rule gatk_recalibrate_primary:
    priority: 1000
    input:
        marked_bam=rules.picard_mark_duplicates.output.marked_bam,
        marked_bai=rules.samtools_index_marked.output,
        genome=ancient("{}.fa".format(project_genome)),
        genome_dict=ancient("{}.dict".format(project_genome)),
        variants=ancient("{}.vcf.gz".format(project_variants)),
        variants_index=ancient("{}.vcf.gz.tbi".format(project_variants))
    output:
        protected("{project_samples}/{sample}/recalibration/{sample}.BaseRecalibrator.primary.grp")
    log:
        "{project_samples}/{sample}/logs/{sample}.BaseRecalibrator.primary.log"
    threads: 1
    resources:
        mem_mb=8192
    shell:
        "gatk --java-options '-Xmx{resources.mem_mb}m' BaseRecalibrator \
        -R {input.genome} \
        -I {input.marked_bam} \
        -known-sites {input.variants} \
        -O {output} 2> {log}"

rule gatk_apply_BQSR:
    priority: 1000
    input:
        marked_bam=rules.picard_mark_duplicates.output.marked_bam,
        marked_bai=rules.samtools_index_marked.output,
        genome=ancient("{}.fa".format(project_genome)),
        genome_dict=ancient("{}.dict".format(project_genome)),
        recal=rules.gatk_recalibrate_primary.output
    output:
        bam=protected("{project_samples}/{sample}/recalibration/{sample}.recalibrated.bam"),
        bai=protected("{project_samples}/{sample}/recalibration/{sample}.recalibrated.bai")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.ApplyBQSR.log")
    params:
        tmp_dir=config["tmp_dir"]
    threads: 1
    resources:
        mem_mb=8192
    shell:
        "gatk  --java-options '-Xmx{resources.mem_mb}m' ApplyBQSR \
        -R {input.genome} \
        -I {input.marked_bam} \
        -bqsr {input.recal} \
        -O {output.bam} \
        --tmp-dir {params.tmp_dir} \
        --static-quantized-quals 10 \
        --static-quantized-quals 15 \
        --static-quantized-quals 20 \
        --static-quantized-quals 25 \
        --static-quantized-quals 30 \
        --static-quantized-quals 35 \
        --static-quantized-quals 40 \
        --static-quantized-quals 45 \
        --static-quantized-quals 50 \
        2> {log}"

rule samtools_index_recal:
    priority: 1000
    input:
        rules.gatk_apply_BQSR.output.bai
    output:
        protected("{project_samples}/{sample}/recalibration/{sample}.recalibrated.bam.bai")
#    params:
#        old_bai="{project_samples}/{sample}/recalibration/{sample}.recalibrated.bai"
    threads: 4
    shell:
        #"samtools index -@ {threads} {input} {output}"
        "mv {input} {output}"

rule gatk_recalibrate_secondary:
    priority: 1000
    input:
        recal_bam=rules.gatk_apply_BQSR.output.bam,
        recal_bai=rules.samtools_index_recal.output,
        genome=ancient("{}.fa".format(project_genome)),
        genome_dict=ancient("{}.dict".format(project_genome)),
        variants=ancient("{}.vcf.gz".format(project_variants)),
        variants_index=ancient("{}.vcf.gz.tbi".format(project_variants)),
        recal_table=rules.gatk_recalibrate_primary.output
    output:
        protected("{project_samples}/{sample}/recalibration/{sample}.BaseRecalibrator.secondary.grp")
    log:
        "{project_samples}/{sample}/logs/{sample}.BaseRecalibrator.secondary.log"
    threads: 1
    resources:
        mem_mb=8192
    shell:
        "gatk --java-options '-Xmx{resources.mem_mb}m' BaseRecalibrator \
        -R {input.genome} \
        -I {input.recal_bam} \
        -known-sites {input.variants} \
        -O {output} 2> {log}"

rule gatk_recalibrate_analyze:
    priority: 1000
    input:
        recal_before=rules.gatk_recalibrate_primary.output,
        recal_after=rules.gatk_recalibrate_secondary.output
    output:
        protected("{project_samples}/{sample}/recalibration/{sample}.BaseRecalibrator.pdf")
    log:
        "{project_samples}/{sample}/logs/{sample}.BaseRecalibrator.AnalyzeCovariates.log"
    threads: 1
    resources:
        mem_mb=8192
    shell:
        "gatk --java-options '-Xmx{resources.mem_mb}m' AnalyzeCovariates \
        -before {input.recal_before} \
        -after {input.recal_after} \
        -plots {output} 2> {log}"

rule gatk_haplotype_caller:
    priority: 900
    input:
        recal_bam=rules.gatk_apply_BQSR.output.bam,
        recal_bai=rules.samtools_index_recal.output,
        genome=ancient("{}.fa".format(project_genome)),
        genome_dict=ancient("{}.dict".format(project_genome)),
        variants=ancient("{}.vcf.gz".format(project_variants)),
        variants_index=ancient("{}.vcf.gz.tbi".format(project_variants)),
    output:
        protected("{project_samples}/{sample}/variants/{sample}.g.vcf")
    log:
        "{project_samples}/{sample}/logs/{sample}.HaplotypeCaller.log"
    params:
        annotation="--annotation-group AS_StandardAnnotation \
         --annotation-group OrientationBiasMixtureModelAnnotation \
         --annotation-group ReducibleAnnotation \
         --annotation-group StandardAnnotation \
         --annotation-group StandardHCAnnotation",
        max_reads=config["gatk_max_reads"],
        report_num_alleles=config["gatk_report_num_alleles"],
        max_num_alleles=config["gatk_max_num_alleles"],
        all_site_pls=config["gatk_all_site_pls"],
        max_assembly_size=config["gatk_max_assembly_size"],
        assembly_padding=config["gatk_assembly_padding"],
        tmp_dir = config["tmp_dir"]
    threads: 4
    resources:
        mem_mb=32768
    shell:
         "gatk --java-options '-Xmx{resources.mem_mb}m' HaplotypeCaller \
         -I {input.recal_bam} \
         -R {input.genome} \
         -O {output} \
         --emit-ref-confidence GVCF \
         --tmp-dir {params.tmp_dir} \
         --dbsnp {input.variants} \
         --max-reads-per-alignment-start {params.max_reads} \
         --annotate-with-num-discovered-alleles {params.report_num_alleles} \
         --max-alternate-alleles {params.max_num_alleles} \
         --max-assembly-region-size {params.max_assembly_size} \
         --assembly-region-padding {params.assembly_padding} \
         --all-site-pls {params.all_site_pls} \
         {params.annotation} \
         2> {log}"

rule gatk_create_gvcf_map_file:
    input:
         expand("{project_samples}/{sample}/variants/{sample}.g.vcf",
                 zip, project_samples=[project_samples, ]*len(samples_names),
                 sample=sorted(samples_names))
    output:
        "{project_variants}/gvcf_map.tsv"
#    log:
#        "{project_samples}/logs/gvcf_map.log"
    run:
        map_file = open("{output}", w)
        for ele in input:
            name = ele.split("/")[-1].split(".")[0]
            print("{}\t{}".format(name, ele), file=map_file)
        map_file.close()

gatk_chroms = " ".join(["-L {} \\" for ele in range(1, 39) ]) + " -L X"

rule gatk_genomic_dbi_import:
    input:
#        "{variants_dir}/gvcf_map.tsv"
        rules.gatk_create_gvcf_map_file.output
    output:
        directory("{variants_dir}/gvcf_db/")
    log:
        "{variants_dir}/GATK.GenomicDbiImport.log"
    params:
#        db_dir="{variants_dir}/gvcf_db",
        intervals=gatk_chroms,
        tmp_dir=config["tmp_dir"]
    threads: threads_max
    resources:
        mem_mb=memory_max
    shell:
       "gatk --java-options '-Xmx{resources.mem_mb}m' GenomicsDBImport \
       --genomicsdb-workspace-path {output} \
       --batch-size {threads} \
       {params.intervals} \
       --sample-name-map {input} \
       --tmp-dir={params.tmp_dir} \
       --reader-threads {threads} \
       2> {logs}"

rule gatk_joint_calling:
    priority: 500
    input:
        dbi=rules.gatk_genomic_dbi_import.output,
        genome=ancient("{}.fa".format(project_genome)),
        genome_dict=ancient("{}.dict".format(project_genome)),
        variants=ancient("{}.vcf.gz".format(project_variants)),
        variants_index=ancient("{}.vcf.gz.tbi".format(project_variants))
    output: "{variants_dir}/var.raw.vcf"
    log: "{variants_dir}/GATK.JointCalling.log"
    params:
        annotation="--annotation-group AS_StandardAnnotation \
         --annotation-group OrientationBiasMixtureModelAnnotation \
         --annotation-group ReducibleAnnotation \
         --annotation-group StandardAnnotation \
         --annotation-group StandardHCAnnotation",
        report_num_alleles=config["gatk_report_num_alleles"],
        max_num_alleles=config["gatk_max_num_alleles"],
        tmp_dir = config["tmp_dir"]
    threads: threads_max
    resources:
        mem_mb=memory_max
    shell:
         "gatk --java-options '-Xmx{resources.mem_mb}m' GenotypeGVCFs \
         -V gendb:/{input} \
         -R {input.genome} \
         -O {output} \
         --tmp-dir {params.tmp_dir} \
         --dbsnp {input.variants} \
         --annotate-with-num-discovered-alleles {params.report_num_alleles} \
         --max-alternate-alleles {params.max_num_alleles} \
         {params.annotation} \
         2> {log}"

rule picard_collect_variant_calling_metrics:
    priority: 500
    input:
        vcf=rules.gatk_joint_calling.output,
        variants=ancient("{}.vcf.gz".format(project_variants)),
        variants_index=ancient("{}.vcf.gz.tbi".format(project_variants))
    output:
        "{variants_dir}/metrics/var.raw.CollectVariantCallingMetrics.txt"
    log:
        protected("{variants_dir}/logs/CollectVariantCallingMetrics.log")
    params:
        picard_jar = picard_jar
    threads: threads_max
    resources:
        mem_mb=memory_max
    shell:
        "java -Xmx{resources.mem_mb}m -jar {params.picard_jar} \
        CollectVariantCallingMetrics \
        -I {input.vcf} \
        --DBSNP {input.variants} \
        -O {output}
        2> {log}"

rule picard_genotype_concordance:
    priority: 500
    input:
        vcf=rules.gatk_joint_calling.output,
        variants=ancient("{}.vcf.gz".format(project_variants)),
        variants_index=ancient("{}.vcf.gz.tbi".format(project_variants))
    output:
        summary="{variants_dir}/metrics/var.raw.vcf.genotype_concordance_summary_metrics",
        detail="{variants_dir}/metrics/var.raw.vcf.genotype_concordance_detail_metrics"
    log: "{variants_dir}/logs/var.raw.Picard_GenotypeConcordance.log"
    params:
        picard_jar·=·picard_jar,
        tmp_dir = config["tmp_dir"],
        prefix="{variants_dir}/metrics/var.raw.vcf"
    threads: threads_max
    resources:
        mem_mb=memory_max
    shell:
         "java·-Xmx{resources.mem_mb}m·-jar·{params.picard_jar} \
         GenotypeConcordance \
         -CV {input.vcf} \
         -TV {input.variants \
         -O {params.prefix} \
         2> {log}"

rule snpeff_annotate:
    priority: 600
    input:
        vcf=rules.gatk_joint_calling.output,
        snpeff_config=snpeff_database + "snpEff.config",
        snpeff_bin=snpeff_database + "/data/" + genome_prefix + "snpEffectPredictor.bin"
    output:
        vcf="{variants_dir}/var.ann.vcf",
        summary="{variants_dir}/metrics/snpeff_summary.csv"
    log:
        "{variants_dir}/logs/var.SnpEff.log"
    params:
        genome_prefix=genome_prefix,
        snpeff_jar=snpeff_jar
    threads: threads_max
    resources:
        mem_mb=memory_max
    shell:
        "$(java·-Xmx{resources.mem_mb}m·-jar·{params.snpeff_jar} \
        -csvStats {output.summary} \
        {params.genome_prefix} {input.vcf} \
        > {output.vcf}) \
        2> {log}"

rule gatk_exclude_strange:
    priority: 500
    input:
        vcf=rules.snpeff_annotate.output.vcf,
        genome=ancient("{}.fa".format(project_genome)),
        genome_dict=ancient("{}.dict".format(project_genome)),
    output:
        "{variants_dir}/var.ann.rest.vcf"
    log:
        "{variants_dir}/logs/var.GATK_SelectSNPs.log"
    threads: 1
    resources: 
        mem_mb=8192
    shell:
        "gatk --java-options '-Xmx{resources.mem_mb}m' SelectVariants \
        -V {input.vcf} \
        -R {input.genome} \
        -O {output} \
        --exclude-non-variants true \
        --select-type-to-exclude SNP \
        --select-type-to-exclude MNP \
        --select-type-to-exclude INDEL \
        2> {log}"

rule gatk_select_snps:
    priority: 500
    input:
        vcf=rules.snpeff_annotate.output.vcf,
        genome=ancient("{}.fa".format(project_genome)),
        genome_dict=ancient("{}.dict".format(project_genome)),
    output:
        temp("{variants_dir}/var.ann.snp.vcf")
    log:
        "{variants_dir}/logs/var.GATK_SelectSNPs.log"
    threads: 1
    resources: 
        mem_mb=8192
    shell:
        "gatk --java-options '-Xmx{resources.mem_mb}m' SelectVariants \
        -V {input.vcf} \
        -R {input.genome} \
        -O {output} \
        --exclude-non-variants true \
        --select-type-to-include SNP \
        --select-type-to-include MNP \
        2> {log}"

rule gatk_filter_snps:
    priority: 500
    input:
        vcf=rules.gatk_select_snps,
        genome=ancient("{}.fa".format(project_genome)),
        genome_dict=ancient("{}.dict".format(project_genome)),
    output:
        protected("{variants_dir}/var.ann.snp.filtered.snp.vcf")
    log:
        "{variants_dir}/logs/var.GATK_FilterSNPs.log"
    threads: 1
    resources:
        mem_mb=8192
    shell:
        "gatk --java-options '-Xmx{resources.mem_mb}m' \
        VariantFiltration \
        -V {input.vcf} \
        -R {input.genome} \
        -O {output} \
        --filterExpression \"\
        QD < 2.0 || \
        FS > 60.0 || \
        MQ < 40.0 || \
        MQRankSum < -12.5 || \
        ReadPosRankSum < -8.0 || \
        SOR > 3.0\" \
        --filterName \"GENERIC_SNP_FILTER\" \
        2> {log}"

rule gatk_select_indels:
    priority: 500
    input:
        vcf=rules.snpeff_annotate.output.vcf,
        genome=ancient("{}.fa".format(project_genome)),
        genome_dict=ancient("{}.dict".format(project_genome)),
    output:
        temp("{variants_dir}/var.ann.indel.vcf")
    log:
        "{variants_dir}/logs/var.GATK_SelectSNPs.log"
    threads: 1
    resources: 
        mem_mb=8192
    shell:
        "gatk --java-options '-Xmx{resources.mem_mb}m' SelectVariants \
        -V {input.vcf} \
        -R {input.genome} \
        -O {output} \
        --exclude-non-variants true \
        --select-type-to-include INDEL \
        2> {log}"

rule gatk_filter_indel:
    priority: 500
    input:
        vcf=rules.gatk_select_indels,
        genome=ancient("{}.fa".format(project_genome)),
        genome_dict=ancient("{}.dict".format(project_genome)),
    output:
        protected("{variants_dir}/var.ann.snp.filtered.snp.vcf")
    log:
        "{variants_dir}/logs/var.GATK_FilterInDels.log"
    threads: 1
    resources:
        mem_mb=8192
    shell:
        "gatk --java-options '-Xmx{resources.mem_mb}m' \
        VariantFiltration \
        -V {input.vcf} \
        -R {input.genome} \
        -O {output} \
        --filterExpression \"\
        QD < 2.0 || \
        FS > 200.0 || \
        ReadPosRankSum < -20.0 || \
        InbreedingCoeff < -0.8 || \
        SOR > 10.0\" \
        --filterName \"GENERIC_INDEL_FILTER\" \
        2> {log}"

rule multiqc_for_intervals:
    priority: 10000
    input:
        hs_metrics=expand("{project_samples}/{sample}/metrics/"
                               "{sample}.HsMetrics.txt",
                               zip,
                               project_samples=[project_samples, ]*len(samples_names),
                               sample=sorted(samples_names)),
    output:
        "{project_main}/MultiQCReportWithIntervals/multiqc_report.html"
    log:
        "{project_main}/logs/multiqc_report_with_intervals.log"
    params:
        input_dir=project_main,
        output_dir="{project_main}/MultiQCReportWithIntervals/".format(project_main=project_main)
    shell:
        "multiqc {params.input_dir} -o {params.output_dir} > {log}"

rule multiqc_for_validation:
    priority: 10000
    input:
        hs_metrics=expand("{project_samples}/{sample}/metrics/"
                               "{sample}.ValidateSamFile.txt",
                               zip,
                               project_samples=[project_samples, ]*len(samples_names),
                               sample=sorted(samples_names)),
    output:
        "{project_main}/MultiQCReportWithValidation/multiqc_report.html"
    log:
        "{project_main}/logs/multiqc_report_with_Validation.log"
    params:
        input_dir=project_main,
        output_dir="{project_main}/MultiQCReportWithValidation/".format(project_main=project_main)
    shell:
        "multiqc {params.input_dir} -o {params.output_dir} > {log}"

rule multiqc:
    priority: 10000
    input:
        gc_bias_metrics=expand("{project_samples}/{sample}/metrics/"
                               "{sample}.GCBiasMetrics.txt",
                               zip,
                               project_samples=[project_samples, ]*len(samples_names),
                               sample=sorted(samples_names)),
        wgs_metrics=expand("{project_samples}/{sample}/metrics/"
                            "{sample}.WgsMetrics.txt",
                            zip,
                            project_samples=[project_samples, ]*len(samples_names),
                            sample=sorted(samples_names)),
        alignment_metrics=expand("{project_samples}/{sample}/metrics/"
                            "{sample}.AlignmentSummaryMetrics.txt",
                            zip,
                            project_samples=[project_samples, ]*len(samples_names),
                            sample=sorted(samples_names)),
        size_metrics=expand("{project_samples}/{sample}/metrics/"
                            "{sample}.InsertSizeMetrics.txt",
                            zip,
                            project_samples=[project_samples, ]*len(samples_names),
                            sample=sorted(samples_names)),
        mark_duplicates=expand("{project_samples}/{sample}/metrics/"
                               "{sample}.MarkDuplicates.txt",
                               zip,
                               project_samples=[project_samples, ]*len(samples_names),
                               sample=sorted(samples_names)),
        idxstats=expand("{project_samples}/{sample}/metrics/{sample}.idxstats",
                        zip,
                        project_samples=[project_samples, ]*len(samples_names),
                        sample=sorted(samples_names)),
        stats=expand("{project_samples}/{sample}/metrics/{sample}.stats",
                     zip,
                     project_samples=[project_samples, ]*len(samples_names),
                     sample=sorted(samples_names)),
        flagstats=expand("{project_samples}/{sample}/metrics/{sample}.flagstats",
                         zip,
                         project_samples=[project_samples, ]*len(samples_names),
                         sample=sorted(samples_names)),
        fastqc_files = expand("{project_samples}/{sample}/metrics/"
                              "fastqc/{sample}_{pair}/{sample}_{pair}_fastqc.zip",
                              zip,
                              project_samples=[project_samples, ]*len(samples_names)*2,
                              sample=sorted(samples_names*2),
                              pair= reads_pairs * len(samples_names)),
        gatk_recal=expand("{project_samples}/{sample}/recalibration/{sample}.BaseRecalibrator.pdf",
                           zip,
                           project_samples=[project_samples, ]*len(samples_names),
                           sample=sorted(samples_names))
    output:
        "{project_main}/MultiQCReport.tar.gz"
    log:
        "{project_main}/logs/multiqc_report.log"
    params:
        input_dir=project_samples,
        output_dir="{project_main}/MultiQCReport/".format(project_main=project_main)
    shell:
        "multiqc {params.input_dir} -o {params.output_dir} 2> {log} && \
        tar -zcvf {output} {params.output_dir} && \
        rm -rf {params.output_dir}"

rule call_variants:
    input:
        #gvcf_map="{variants_dir}/gvcf_map.tsv".format(variants_dir=variants_dir),
        var_filtered_snp="{variants_dir}/var.ann.snp.filtered.vcf",
        var_filtered_indel="{variants_dir}/var.ann.indel.filtered.vcf",
        var_rest="{variants_dir}/var.ann.rest.vcf",
        picard_genotype_concordance_summary="{variants_dir}/metrics/var.raw.vcf.genotype_concordance_summary_metrics",
        picard_genotype_concordance_detail="{variants_dir}/metrics/var.raw.vcf.genotype_concordance_detail_metrics",
        picard_variant_calling_metrics="{variants_dir}/metrics/var.raw.CollectVariantCallingMetrics.txt",
        #var_ann="{variants_dir}/var.ann.vcf".format(variants_dir=variants_dir),
        #picard_genotype_concordance_summary="{variants_dir}/metrics/var.raw.vcf.genotype_concordance_summary_metrics".format(variants_dir=variants_dir),
        #picard_genotype_concordance_detail="{variants_dir}/metrics/var.raw.vcf.genotype_concordance_detail_metrics".format(variants_dir=variants_dir),
        #picard_variant_calling_metrics="{variants_dir}/metrics/var.raw.CollectVariantCallingMetrics.txt".format(variants_dir=variants_dir),
    output:
        directory("{variants_dir}/MultiQCReport/")
    log:
        "{variants_dir}/logs/multiqc_report.log"
    params:
        input_dir=variants_dir
    shell:
        "multiqc {params.input_dir} -o {output} 2> {log}"

rule compress_multiqc_variant:
    priority: 20000
    input:
        rules.call_variants.output
    output:
        "{variants_dir}/MultiQCReport.tar.gz".format(variants_dir=variants_dir)
    shell:
        "tar·-zcvf·{output}·{input}"

