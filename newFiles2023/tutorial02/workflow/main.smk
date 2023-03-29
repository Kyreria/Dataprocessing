SAMPLES = ["A", "B", "C"]

rule all:
    input:
        "results/out.html"

rule bwa_map:
    input:
        "resources/genome.fa",
        "resources/samples/{sample}.fastq"
    output:
        "results/mapped_reads/{sample}.bam"
    message: "Executing bwa mem on the following {input} to generate the following {output}"
    shell:
        "bwa mem {input} | samtools view -Sb -> {output}"

rule samtools_sort:
    input:
        "results/mapped_reads/{sample}.bam"
    output:
        "results/sorted_reads/{sample}.bam"
    message:"Executing samtools on the following {input} to generate the following {output}"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"

rule samtools_index:
    input:
        "results/sorted_reads/{sample}.bam"
    output:
        "results/sorted_reads/{sample}.bam.bai"
    message: "Indexing {input} with samtools to the following {output}"
    shell:
        "samtools index {input}"

rule bcftools_call:
    input:
        fa="resources/genome.fa",
        bam=expand("results/sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("results/sorted_reads/{sample}.bam", sample=SAMPLES)
    output:
        "results/calls/all.vcf"
    message: ""
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv -> {output}"

rule report:
    input:
        "results/calls/all.vcf"
    output:
        "results/out.html"
    run:
        from snakemake.utils import report
        with open(input[0]) as f:
            n_calls = sum(1 for line in f if not line.startswith("#"))

        report("""
        An example workflow
        ==================================
        
        Reads were mapped to the Yeas reference genome
        and variants were called jointly with
        SAMtools/BCFtools.
        
        This resulted in {n_calls} variants (see Table T1_).
        """, output[0], metadata="Author : Dennis Haandrikman", T1=input[0])