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

rule all:
    input:
        multiqc="{project_main}/MultiQCReport/multiqc_report.html".format(project_main=project_main),
        gatk_recal=expand("{project_samples}/{sample}/recalibration/{sample}.recal.bam.bai",
                          zip,
                          project_samples=[project_samples, ]*len(samples_names),
                          sample=sorted(samples_names))

rule samtools_index_recal:
    input:
        "{project_samples}/{sample}/recalibration/{sample}.recal.bam"
    output:
        protected("{project_samples}/{sample}/recalibration/{sample}.recal.bam.bai")
    threads: 20
    shell:
        "samtools index -@ {threads} {input} {output}"

rule gatk_apply_BQSR:
    input:
        marked_bam="{project_samples}/{sample}/{sample}.bam",
        marked_bai="{project_samples}/{sample}/{sample}.bai",
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
        "gatk  --java-options {params.memory} ApplyBQSR \
        -R {input.genome} \
        -I {input.marked_bam} \
        -bqsr {input.recal} \
        -O {output} \
        -SQQ 10 -SQQ 15 -SQQ 20 -SQQ 25 -SQQ 30 -SQQ 35 -SQQ 40 -SQQ 45 -SQQ 50"

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
        "gatk --java-options {params.memory} AnalyzeCovariates \
        -before {input.before} \
        -after {input.after} \
        -plots {output} > {log}"

rule gatk_recalibrate_2nd:
    input:
        marked_bam="{project_samples}/{sample}/{sample}.bam",
        marked_bai="{project_samples}/{sample}/{sample}.bai",
        genome=ancient("{}.fa".format(project_genome)),
        genome_dict=ancient("{}.dict".format(project_genome)),
        variants=ancient("{}.vcf.gz".format(project_variants)),
        recal_table="{project_samples}/{sample}/recalibration/{sample}.recal.1st.table"
    output:
        protected("{project_samples}/{sample}/recalibration/{sample}.recal.2nd.table")
    log:
        "{project_samples}/{sample}/logs/{sample}.recal.2nd.log"
    threads: 20
    params:
        memory="-Xmx160g"
    shell:
        "gatk --java-options {params.memory} BaseRecalibrator \
        -R {input.genome} \
        -I {input.marked.bam} \
        -knownSites {input.variants} \
        –bqsr {input.recal_table} \
        -O {output} > {log}"

rule gatk_recalibrate_1st:
    input:
        marked_bam="{project_samples}/{sample}/{sample}.bam",
        marked_bai="{project_samples}/{sample}/{sample}.bai",
        genome=ancient("{}.fa".format(project_genome)),
        genome_dict=ancient("{}.dict".format(project_genome)),
        variants=ancient("{}.vcf.gz".format(project_variants))
    output:
        protected("{project_samples}/{sample}/recalibration/{sample}.recal.1st.table")
    log:
        "{project_samples}/{sample}/logs/{sample}.recal.1st.log"
    threads: 20
    params:
        memory="-Xmx160g"
    shell:
        "gatk --java-options {params.memory} BaseRecalibrator \
        -R {input.genome} \
        -I {input.marked.bam} \
        -knownSites {input.variants} \
        -O {output} > {log}"

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
        gatk_recal=expand("{project_samples}/{sample}/recalibration/{sample}.recal.pdf",
                          zip,
                          project_samples=[project_samples, ]*len(samples_names),
                          sample=sorted(samples_names))
    output:
        "{project_main}/MultiQCReport/multiqc_report.html"
    log:
        "{project_main}/logs/multiqc_report.log"
    params:
        output_dir="{project_main}/MultiQCReport/".format(project_main=project_main)
    shell:
        "multiqc {project_samples} -o {params.output_dir} > {log}"

rule picard_validate_sam_file:
    input:
        marked_bam="{project_samples}/{sample}/{sample}.bam",
        marked_bai="{project_samples}/{sample}/{sample}.bai",
        genome="{}.fa".format(project_genome),
        genome_dict="{}.dict".format(project_genome)
    output:
        txt=protected("{project_samples}/{sample}/metrics/{sample}.ValidateSamFile.txt")
    log:
        protected("{project_samples}/{sample}/metrics/{sample}.ValidateSamFile.log")
    threads: 4
    params:
        picard_jar = picard_jar
    shell:
        "java -Xmx32g -jar {params.picard_jar} ValidateSamFile \
        I={input.marked_bam} OUTPUT={output.txt} \
        R={input.genome} MODE='VERBOSE' > {log}"

rule picard_gc_bias_metrics:
    input:
        marked_bam="{project_samples}/{sample}/{sample}.bam",
        marked_bai="{project_samples}/{sample}/{sample}.bam.bai",
        genome="{}.fa".format(project_genome),
        genome_dict="{}.dict".format(project_genome)
    output:
        metrics=protected("{project_samples}/{sample}/metrics/{sample}.GCBiasMetrics.txt"),
        chart=protected("{project_samples}/{sample}/metrics/{sample}.GCBiasMetrics.pdf"),
        summary=protected("{project_samples}/{sample}/metrics/{sample}.GCBiasMetricsSummary.txt")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.GCBiasMetrics.log")
    threads: 4
    params:
        picard_jar = picard_jar
    shell:
        "java -Xmx32g -jar {params.picard_jar} CollectGcBiasMetrics \
        R={input.genome} I={input.marked_bam} O={output.metrics} S={output.summary} \
        CHART={output.chart} > {log}"

rule picard_wgs_metrics:
    input:
        marked_bam="{project_samples}/{sample}/{sample}.bam",
        marked_bai="{project_samples}/{sample}/{sample}.bam.bai",
        genome="{}.fa".format(project_genome),
        genome_dict="{}.dict".format(project_genome)
    output:
        txt=protected("{project_samples}/{sample}/metrics/{sample}.WgsMetrics.txt")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.WgsMetrics.log")
    threads: 4
    params:
        picard_jar = picard_jar
    shell:
        "java -Xmx32g -jar {params.picard_jar} CollectWgsMetrics \
        R={input.genome} I={input.marked_bam} O={output.txt} > {log}"

rule picard_alignment_summary:
    input:
        marked_bam="{project_samples}/{sample}/{sample}.bam",
        marked_bai="{project_samples}/{sample}/{sample}.bam.bai",
        genome="{}.fa".format(project_genome),
        genome_dict="{}.dict".format(project_genome)
    output:
        txt=protected("{project_samples}/{sample}/metrics/{sample}.AlignmentSummaryMetrics.txt")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.AlignmentSummaryMetrics.log")
    threads: 4
    params:
        picard_jar = picard_jar
    shell:
        "java -Xmx32g -jar {params.picard_jar} CollectAlignmentSummaryMetrics \
        R={input.genome} I={input.marked_bam} O={output.txt} > {log}"

rule picard_size_metrics:
    input:
        marked_bam="{project_samples}/{sample}/{sample}.bam",
        marked_bai="{project_samples}/{sample}/{sample}.bam.bai",
    output:
        txt=protected("{project_samples}/{sample}/metrics/{sample}.InsertSizeMetrics.txt"),
        pdf=protected("{project_samples}/{sample}/metrics/{sample}.InsertSizeMetrics.pdf")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.InsertSizeMetrics.log")
    params:
        picard_jar = picard_jar
    threads: 4
    shell:
        "java -Xmx32g -jar {params.picard_jar} CollectInsertSizeMetrics \
        I={input.marked_bam} \
        O={output.txt} H={output.pdf} M=0.5 > {log}"

rule samtools_stats:
    input:
        marked_bam="{project_samples}/{sample}/{sample}.bam",
        marked_bai="{project_samples}/{sample}/{sample}.bam.bai"
    output:
        protected("{project_samples}/{sample}/metrics/{sample}.stats")
    threads: 1
    shell:
        "samtools stats {input.marked_bam} > {output}"

rule samtools_idxstats:
    input:
        marked_bam="{project_samples}/{sample}/{sample}.bam",
        marked_bai="{project_samples}/{sample}/{sample}.bam.bai"
    output:
        protected("{project_samples}/{sample}/metrics/{sample}.idxstats")
    threads: 1
    shell:
        "samtools idxstats {input.marked_bam} > {output}"

rule samtools_flagstats:
    input:
        marked_bam="{project_samples}/{sample}/{sample}.bam",
        marked_bai="{project_samples}/{sample}/{sample}.bam.bai"
    output:
        protected("{project_samples}/{sample}/metrics/{sample}.flagstats")
    threads: 20
    shell:
        "samtools flagstat -@ {threads} {input.marked_bam} > {output}"

rule samtools_index_marked:
    input:
        "{project_samples}/{sample}/{sample}.bam"
    output:
        protected("{project_samples}/{sample}/{sample}.bam.bai")
    threads: 20
    shell:
        "samtools index -@ {threads} {input} {output}"

rule picard_mark_duplicates:
    input:
        sorted_bam="{project_samples}/{sample}/{sample}.sorted.bam",
        sorted_bai="{project_samples}/{sample}/{sample}.sorted.bam.bai"
    output:
        marked_bam=protected("{project_samples}/{sample}/{sample}.bam"),
        marked_metrics=protected("{project_samples}/{sample}/metrics/{sample}.MarkDuplicates.txt")
    log:
        protected("{project_samples}/{sample}/logs/{sample}.MarkDuplicates.log")
    params:
        picard_jar=picard_jar
    threads: 20
    shell:
        "java -Xmx160g -jar {params.picard_jar} \
         MarkDuplicates I={input.sorted_bam} \
         O={output.marked_bam} \
         M={output.marked_metrics} > {log}"

rule samtools_index_sorted:
    input:
        "{project_samples}/{sample}/{sample}.sorted.bam"
    output:
        temp("{project_samples}/{sample}/{sample}.sorted.bam.bai")
    threads: 20
    shell:
        "samtools index -@ {threads} {input} {output}"

rule samtools_sort:
    input:
        "{project_samples}/{sample}/{sample}.unsorted.bam"
    output:
        temp("{project_samples}/{sample}/{sample}.sorted.bam")
    threads: 20
    params:
        memory="4G"
    shell:
        "samtools sort {input} -o {output} -@ {threads} -m {params.memory}"

# rule picard_add_metadata:
#     input:
#         "{project_samples}/{sample}/{sample}.missing.bam"
#     output:
#         temp("{project_samples}/{sample}/{sample}.unsorted.bam")
#     params:
#         picard_jar = picard_jar,
#         rgsm = "{sample}"
#     threads: 4
#     shell:
#         "java -Xmx32g -jar {params.picard_jar} AddOrReplaceReadGroups \
#         I={input} O={output} \
#         RGID={params.rgsm} \
#         RGSM={params.rgsm} \
#         RGPL=illumina \
#         RGLB={params.rgsm} \
#         RGPU={params.rgsm}"

rule bwa_map_reads:
    input:
        reads_1=ancient("{project_samples}/{sample}/{sample}_1.fq.gz"),
        reads_2=ancient("{project_samples}/{sample}/{sample}_2.fq.gz"),
        genome="{project_genome}.ann".format(project_genome=project_genome)
    output:
        # temp("{project_samples}/{sample}/{sample}.missing.bam")
        temp("{project_samples}/{sample}/{sample}.unsorted.bam")
    params:
        prefix="{project_genome}".format(project_genome=project_genome),
        rg=r"@RG\tID:{sample}\tSM:{sample}\tPL:illumina\tLB:{sample}\tPU:{sample}"
    threads:20
    log:
        "{project_samples}/{sample}/{sample}.bwa.log"
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {params.prefix} \
        {input.reads_1} {input.reads_2} | samtools view -Sb - > {output}) 2> {log}"

rule prepare_genome_index:
    input:
        ancient("{project_genome}.fa")
    output:
        faidx="{project_genome}.fa.fai",
        dict="{project_genome}.dict"
    params:
        picard_jar = picard_jar
    run:
        shell("samtools faidx {input} -o {output.faidx}")
        shell("java -jar {params.picard_jar} "
              "CreateSequenceDictionary "
              "R={input} "
              "O={output.dict}")

rule prepare_bwa_genome:
    input:
        ancient("{project_genome}.fa")
    output:
        "{project_genome}.ann"
    params:
        picard_jar = picard_jar
    shell:
        "bwa index -p {project_genome} {input}"

rule fastqc:
    input:
        ancient("{project_samples}/{sample}/{sample}_{pair}.fq.gz")
    output:
        "{project_samples}/{sample}/metrics/fastqc/{sample}_{pair}/{sample}_{pair}_fastqc.zip",
    params:
        output_dir="{project_samples}/{sample}/metrics/fastqc/{sample}_{pair}"
    shell:
        "fastqc {input} -o {params.output_dir}"
