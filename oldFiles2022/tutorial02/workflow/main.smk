SAMPLES = []

rule all:
    input:
        'resources/report/out.html'

rule bwa_map:
    input:
        "resources/genome.fa",
        "resources/samples/{sample}.fastq"
    output:
        "resources/mapped_reads/{sample}.bam"
    message: "Executing bwa mem on the following {input} to generate the following {output}"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"

rule samtools_sort:
    input:
        "resources/mapped_reads/{sample}.bam"
    output:
        "resources/sorted_reads/{sample}.bam"
    message: "Sorting the reads with of the following {input} to generate the following sorted {output}"
    shell:
        "samtools sort -T resources/sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"

rule samtools_index:
    input:
        "resources/sorted_reads/{sample}.bam"
    output:
        "resources/sorted_reads/{sample}.bam.bai"
    message: "Indexing the sorted reads of {input} into the following {output}"
    shell:
        "samtools index {input}"

rule bcftools_call:
    input:
        fa="resources/genome.fa",
        bam=expand("resources/sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("resources/sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "resources/calls/all.vcf"
    message: "Calling out variants of the mapped reads into {output}"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"

rule report:
    input:
        "resources/calls/all.vcf"
    output:
        "resources/report/out.html"
    message: "Writing a report with the following {input} into {output}"
    run:
        from snakemake.utils import report
        with open(input[0]) as f:
            n_calls = sum(1 for line in f if not line.startswith("#"))

        report("""
        An example workflow
        ===================================
        
        Reads were mapped to the Yeas reference genome 
        and variants were called jointly with
        SAMtools/BCFtools.
        
        This resulted in {n_calls} variants (see Table T1_).
        """, output[0], metadata="Author: Mr Pipeline", T1=input[0])

rule dag_file:
    input:
        "resources/genome.fa"