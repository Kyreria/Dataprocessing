SAMPLES = ['Sample1', 'Sample2', 'Sample3']

rule all:
    input:
        'results/merged.txt'

rule quantify_genes:
    input:
        genome = 'resources/genome.fa',
        r1 = 'resources/{sample}.R1.fastq.gz',
        r2 = 'resources/{sample}.R2.fastq.gz'
    output:
        temp('results/{sample}.txt')
    shell:
        'echo {input.genome} {input.r1} {input.r2} > {output}'

rule merged_results:
    input:
        expand('results/{sample}.txt',sample=SAMPLES)
    output:
        'results/merged.txt'
    shell:
        'cat {input} > {output}'

rule clean:
    shell:
        'rm results/*.txt'