# Snakefile for mapping reads from Whole Genome NGS

# All input parameters are being read from config.yaml
# Script will work on files in "project_dir"
# genome in FASTA format is located in genome_dir and is prefixed with
# "genome_prefix"like:
# It assumes reads are located in "reads_dir" and they have the structure:
# {sample}/{file_name}_{pair}.fq.gz, where {pair} may be "1" or "2"

# Selected Picard command to implement
# java -Xmx32g -jar /opt/conda/share/picard-2.18.16-0/picard.jar CollectHsMetrics I=CFA_615.bam O=CFA_615.hs_metrics.txt R=/data/OPUS/genome/CanFam3.1_Ensembl94.fa BAIT_INTERVALS=/data/OPUS/genome/CanFam3.1_Ensembl94.intervals TARGET_INTERVALS=/data/OPUS/genome/CanFam3.1_Ensembl94.intervals

project_main = config["project_dir"]
project_reads = "/".join((project_main, config["reads_dir"]))
project_genome = "/".join((project_main, config["genome_dir"], config["genome_prefix"]))
project_variants = "/".join((project_main, config["variants_dir"], config["genome_prefix"]))
project_samples = "/".join((project_main, config["samples_dir"]))
reads_pairs = ["1", "2"]

(picard_version) = glob_wildcards("/opt/conda/share/picard-{version}/picard.jar")
if not len(picard_version.version):
    raise ValueError("no Picard jar file found")
picard_jar = "/opt/conda/share/picard-{}/picard.jar".format(picard_version.version[0])

print("Using Picard version {}".format(picard_jar))

def reads_files_group():
    (dirs, ) = glob_wildcards(project_samples + "/{dir}")
    dirs = [ele for ele in dirs if "/" not in ele]
    return dirs

samples_names = reads_files_group()

print("{} samples loaded".format(len(samples_names)))
print(*samples_names, sep=", ")

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
        ancient("{project_genome}.fa")
    output:
        protected("{project_genome}.fa.fai")
    run:
        "samtools faidx {input} -o {output}"

rule prepare_genome_dictionary:
    input:
        fa=ancient("{project_genome}.fa"),
        fai=rules.prepare_genome_index.output
    output:
        protected("{project_genome}.dict")
    params:
        picard_jar = picard_jar
    threads: 1
    run:
        "java -jar {params.picard_jar} CreateSequenceDictionary \
        R={input} O={output.dict}"

rule prepare_bwa_genome:
    input:
        ancient("{project_genome}.fa")
    output:
        amb=protected("{project_genome}.amb"),
        ann=protected("{project_genome}.ann"),
        bwt=protected("{project_genome}.bwt"),
        pac=protected("{project_genome}.pac"),
        sa=protected("{project_genome}.sa")
    threads: 1
    shell:
        "bwa index -p {project_genome} {input}"

rule bwa_map_reads:
    priority: 1000
    input:
        reads_1=ancient("{project_samples}/{sample}/{sample}_1.fq.gz"),
        reads_2=ancient("{project_samples}/{sample}/{sample}_2.fq.gz"),
        genome=rules.prepare_bwa_genome.output
    output:
        temp("{project_samples}/{sample}/{sample}.unsorted.bam")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.bwa.log")
    params:
        prefix="{project_genome}".format(project_genome=project_genome),
        rg=r"@RG\tID:{sample}\tSM:{sample}\tPL:illumina\tLB:{sample}\tPU:{sample}"
    threads: threads_max
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {params.prefix} \
        {input.reads_1} {input.reads_2} | samtools view -Sb - > {output}) 2> {log}"

rule samtools_sort_bwa_map:
    priority: 1000
    input:
        rules.bwa_map_reads.output
    output:
        temp("{project_samples}/{sample}/{sample}.sorted.bam")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.sort_bwa_map.log")
    threads: threads_max
    resources:
        mem_mb="4096M"
    shell:
        "samtools sort {input} \
        -o {output} \
        -@ {threads} \
        -m {resources.mem_mb} \
        -T /tmp/{wildcards.sample} 2> {log}"

rule samtools_index_sorted_bwa_map:
    priority: 1000
    input:
        rules.samtools_sort_bwa_map.output
    output:
        temp("{project_samples}/{sample}/{sample}.sorted.bam.bai")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.index_sorted_bwa_map.log")
    threads: threads_max
    shell:
        "samtools index -@ {threads} {input} {output}"

rule picard_mark_duplicates:
    priority: 900
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
        mem_mb = "196608M"
    shell:
        "java -Xmx{resources.mem_mb} -jar {params.picard_jar} \
         MarkDuplicates I={input.sorted_bam} \
         O={output.marked_bam} \
         M={output.marked_metrics} 2> {log}"

rule samtools_index_marked:
    priority: 900
    input:
        rules.picard_mark_duplicates.output.marked_bam
    output:
        protected("{project_samples}/{sample}/{sample}.bam.bai")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.index_marked.log")
    threads: threads_max
    shell:
        "samtools index -@ {threads} {input} {output}"

rule samtools_flagstats:
    priority: 750
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
    priority: 750
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
    priority: 750
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
    priority: 500
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
    threads: 2
    resources:
        mem_mb = 16384
    shell:
        "java -Xmx{resources.mem_mb} -jar {params.picard_jar} \
        CollectInsertSizeMetrics \
        I={input.marked_bam} \
        O={output.txt} \
        H={output.pdf} \
        M=0.5 2> {log}"

rule picard_alignment_summary:
    priority: 500
    input:
        marked_bam=rules.picard_mark_duplicates.output.marked_bam,
        marked_bai=rules.samtools_index_marked.output
        genome=ancient("{project_genome}.fa"),
        genome_dict=rules.prepare_genome_dictionary.output
    output:
        txt=protected("{project_samples}/{sample}/metrics/{sample}.AlignmentSummaryMetrics.txt")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.AlignmentSummaryMetrics.log")
    params:
        picard_jar = picard_jar
    threads: 2
    resources:
        mem_mb = 16384
    shell:
        "java -Xmx{resources.mem_mb} -jar {params.picard_jar} \
        CollectAlignmentSummaryMetrics \
        R={input.genome} \
        I={input.marked_bam} \
        O={output.txt} 2> {log}"

rule picard_wgs_metrics:
    priority: 500
    input:
        marked_bam=rules.picard_mark_duplicates.output.marked_bam,
        marked_bai=rules.samtools_index_marked.output
        genome=ancient("{project_genome}.fa"),
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
        "java -Xmx{resources.mem_mb} -jar {params.picard_jar} \
        CollectWgsMetrics \
        R={input.genome} \
        I={input.marked_bam} \
        O={output.txt} 2> {log}"

rule picard_gc_bias_metrics:
    priority: 500
    input:
        marked_bam=rules.picard_mark_duplicates.output.marked_bam,
        marked_bai=rules.samtools_index_marked.output
        genome=ancient("{project_genome}.fa"),
        genome_dict=rules.prepare_genome_dictionary.output
    output:
        metrics=protected("{project_samples}/{sample}/metrics/{sample}.GCBiasMetrics.txt"),
        chart=protected("{project_samples}/{sample}/metrics/{sample}.GCBiasMetrics.pdf"),
        summary=protected("{project_samples}/{sample}/metrics/{sample}.GCBiasMetricsSummary.txt")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.GCBiasMetrics.log")
    params:
        picard_jar = picard_jar
    threads: 2
    resources:
        mem_mb = 16384
    shell:
        "java -Xmx{resources.mem_mb} -jar {params.picard_jar} \
        CollectGcBiasMetrics \
        R={input.genome} I={input.marked_bam} O={output.metrics} S={output.summary} \
        CHART={output.chart} 2> {log}"

rule picard_validate_sam_file:
    priority: 500
    input:
        marked_bam=rules.picard_mark_duplicates.output.marked_bam,
        marked_bai=rules.samtools_index_marked.output
        genome=ancient("{project_genome}.fa"),
        genome_dict=rules.prepare_genome_dictionary.output
    output:
        txt=protected("{project_samples}/{sample}/metrics/{sample}.ValidateSamFile.txt")
    log:
        protected("{project_samples}/{sample}/metrics/{sample}.ValidateSamFile.log")
    params:
        picard_jar = picard_jar
    threads: 2
    resources:
        mem_mb = 16384
    shell:
        "java -Xmx{resources.mem_mb} -jar {params.picard_jar} \
        ValidateSamFile \
        I={input.marked_bam} OUTPUT={output.txt} \
        R={input.genome} MODE='VERBOSE' 2> {log}"

rule multiqc:
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
        # gatk_recal=expand("{project_samples}/{sample}/recalibration/{sample}.recal.pdf",
        #                   zip,
        #                   project_samples=[project_samples, ]*len(samples_names),
        #                   sample=sorted(samples_names))
    output:
        "{project_main}/MultiQCReport/multiqc_report.html"
    log:
        "{project_main}/logs/multiqc_report.log"
    params:
        output_dir="{project_main}/MultiQCReport/".format(project_main=project_main)
    shell:
        "multiqc {project_samples} -o {params.output_dir} > {log}"

rule gatk_index_variants:
    input:
        "{}.vcf.gz".format(project_variants)
    output:
        "{}.vcf.gz.tbi".format(project_variants)
    log:
        "{}.log".format(project_variants)
    params:
        memory="-Xmx160g"
    shell:
        "gatk IndexFeatureFile \
        -F {input} \
        -O {output} 2> {log}"

rule gatk_recalibrate_1st:
    input:
        marked_bam="{project_samples}/{sample}/{sample}.bam",
        marked_bai="{project_samples}/{sample}/{sample}.bam.bai",
        genome=ancient("{}.fa".format(project_genome)),
        genome_dict=ancient("{}.dict".format(project_genome)),
        variants=ancient("{}.vcf.gz".format(project_variants)),
        variants_index=ancient("{}.vcf.gz.tbi".format(project_variants))
    output:
        protected("{project_samples}/{sample}/recalibration/{sample}.recal.1st.table")
    log:
        "{project_samples}/{sample}/logs/{sample}.recal.1st.log"
    threads: 20
    params:
        memory="-Xmx160g"
    shell:
        "gatk --java-options '{params.memory}' BaseRecalibrator \
        -R {input.genome} \
        -I {input.marked_bam} \
        -known-sites {input.variants} \
        -O {output} 2> {log}"

rule gatk_recalibrate_2nd:
    input:
        marked_bam="{project_samples}/{sample}/{sample}.bam",
        marked_bai="{project_samples}/{sample}/{sample}.bam.bai",
        genome=ancient("{}.fa".format(project_genome)),
        genome_dict=ancient("{}.dict".format(project_genome)),
        variants=ancient("{}.vcf.gz".format(project_variants)),
        variants_index=ancient("{}.vcf.gz.tbi".format(project_variants)),
        recal_table="{project_samples}/{sample}/recalibration/{sample}.recal.1st.table"
    output:
        protected("{project_samples}/{sample}/recalibration/{sample}.recal.2nd.table")
    log:
        "{project_samples}/{sample}/logs/{sample}.recal.2nd.log"
    threads: 20
    params:
        memory="-Xmx160g"
    shell:
        "gatk --java-options '{params.memory}' BaseRecalibrator \
        -R {input.genome} \
        -I {input.marked_bam} \
        -known-sites {input.variants} \
        –BQSR {input.recal_table} \
        -O {output} 2> {log}"

rule gatk_recalibrate_analyze:
    input:
        recal_before="{project_samples}/{sample}/recalibration/{sample}.recal.1st.table",
        recal_after="{project_samples}/{sample}/recalibration/{sample}.recal.2nd.table"
    output:
        protected("{project_samples}/{sample}/recalibration/{sample}.recal.pdf")
    threads: 20
    params:
        memory="-Xmx160g"
    shell:
        "gatk --java-options '{params.memory}' AnalyzeCovariates \
        -before {input.recal_before} \
        -after {input.recal_after} \
        -plots {output} 2> {log}"

rule gatk_apply_BQSR:
    input:
        marked_bam="{project_samples}/{sample}/{sample}.bam",
        marked_bai="{project_samples}/{sample}/{sample}.bam.bai",
        genome=ancient("{}.fa".format(project_genome)),
        genome_dict=ancient("{}.dict".format(project_genome)),
        recal="{project_samples}/{sample}/recalibration/{sample}.recal.1st.table"
    output:
        protected("{project_samples}/{sample}/recalibration/{sample}.recal.bam")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.recal.log")
    threads: 20
    params:
        memory="-Xmx160g"
    shell:
        "gatk  --java-options '{params.memory}' ApplyBQSR \
        -R {input.genome} \
        -I {input.marked_bam} \
        -bqsr {input.recal} \
        -O {output} \
        -SQQ 10 -SQQ 15 \
        -SQQ 20 -SQQ 25 \
        -SQQ 30 -SQQ 35 \
        -SQQ 40 -SQQ 45 \
        -SQQ 50 2> {log}"

rule samtools_index_recal:
    input:
        "{project_samples}/{sample}/recalibration/{sample}.recal.bam"
    output:
        protected("{project_samples}/{sample}/recalibration/{sample}.recal.bam.bai")
    threads: 20
    shell:
        "samtools index -@ {threads} {input} {output}"

rule all:
    input:
        multiqc="{project_main}/MultiQCReport/multiqc_report.html".format(project_main=project_main),
        gatk_recal=expand("{project_samples}/{sample}/recalibration/{sample}.recal.bam.bai",
                          zip,
                          project_samples=[project_samples, ]*len(samples_names),
                          sample=sorted(samples_names))
