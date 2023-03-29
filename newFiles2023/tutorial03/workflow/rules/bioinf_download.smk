url = "bioinf.nl/~fennaf/snakemake/WC03/data/"
# Website seems to be down?
# Skipping this part of the assignment, need to discuss.
# Keep getting 404 errors regarding the page from the HTTP requests.

rule download:
    output:
       "resources/bioinf/{sample}.bam"
    message: "Downloading data from the following link {url}"
    shell:
        "wget {url}{output}"