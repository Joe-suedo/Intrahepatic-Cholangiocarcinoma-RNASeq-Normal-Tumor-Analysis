# 1. Define your dictionary for renaming the files grouped by sample accessions
title_to_samn = {
    '090608N':'SAMN09941351','090608T':'SAMN09941350','098699N':'SAMN09941349',
    '098699T':'SAMN09941348','099535N':'SAMN09941347','099535T':'SAMN09941343',
    '11770N':'SAMN09941342','11770T':'SAMN09941341','81023N':'SAMN09941340',
    '81023T':'SAMN09941346','81032N':'SAMN09941345','81032T':'SAMN09941344',
    '85421N':'SAMN09941337','85421T':'SAMN09941336','86279N':'SAMN09941335',
    '86279T':'SAMN09941334','86960N':'SAMN09941339','86960T':'SAMN09941338',
    '87823N':'SAMN09941333','87823T':'SAMN09941332','88659N':'SAMN09941331',
    '88659T':'SAMN09941330','88674N':'SAMN09941329','88674T':'SAMN09941328',
    '88724N':'SAMN09941327','88724T':'SAMN09941326','89290N':'SAMN09941325',
    '89290T':'SAMN09941324','96484N':'SAMN09941323','96484T':'SAMN09941322'
}

# 2. Rule All (The target)
rule all:
    input:
        "fastqc_1.html"

# 3. Download Rule
rule download_raw_fastq:
    output:"fastq-raw/.download_finished" #Use a file instead of a directory
    conda: "snakemake_env"
    shell:
        "fastq-dl --accession PRJNA488803 --cpus 40 --group-by-sample --outdir fastq-raw"

# 4. Rename Rule
rule rename_raw_fastq:
    input:
        # This MUST match the output of the download rule exactly
        "fastq-raw/fastq-run-mergers.tsv"
    output:
        "fastq-raw/{title}_{read}.fastq.gz"
    run:
        # 1. Get the old SAMN ID from your dictionary
        old_id = title_to_samn[wildcards.title]
        
        # 2. Construct the path to the file fastq-dl created
        # Check if fastq-dl adds _1 or _R1 (most add _1 or _2)
        old_path = f"fastq-raw/{old_id}_{wildcards.read}.fastq.gz"
        
        # 3. Move it to the new name
        shell(f"mv {old_path} {output}")


# 5. FastQC Rule
rule fastqc_before:
    input:
        "fastq-raw/{sample}_{read}.fastq.gz"
    output:
        html="fastqc_1/{sample}_{read}_fastqc.html",
        zip="fastqc_1/{sample}_{read}_fastqc.zip"
    threads: 10
    conda: "qc_env"
    shell:
        "fastqc -t {threads} -o fastqc_1/ {input}"

# 6. MultiQC Rule
rule multiqc:
    input:
        expand("fastqc_1/{title}_{read}_fastqc.html", title=title_to_samn.keys(), read=["R1", "R2"])
    output:
        "fastqc_1.html"
    conda: "qc_env"
    shell:
        "multiqc fastqc_1/ -n fastqc_1 -o ./ "
