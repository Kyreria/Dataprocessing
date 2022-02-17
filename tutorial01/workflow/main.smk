SAMPLES = ['Sample1', 'Sample2', 'Sample3']

rule all:
    input:
        expand('{sample}.txt', sample = SAMPLES)

rule quantify_genes:
    input:
        genome = 'resources/genome.fa',
        r1 = 'resources/{sample}.R1.fastq.gz',
        r2 = 'resources/{sample}.R2.fastq.gz'
    output:
        'results/{sample}.txt'
    shell:
        'echo {input.genome} {input.r1} {input.r2} > {output}'
