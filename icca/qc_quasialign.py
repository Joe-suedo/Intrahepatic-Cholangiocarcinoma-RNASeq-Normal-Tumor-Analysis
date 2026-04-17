# 1. Automatically find samples in fastq-raw directory
SAMPLES, = glob_wildcards("fastq-raw/{sample}_R1.fastq.gz")

# 2. Target Rule: Tells Snakemake exactly what files we want at the end
rule all:
    input:
        expand("quant/{sample}/quant.sf", sample=SAMPLES),
        "fastqc_2.html"

# 3. Download & Index rRNA (Fixed: now inside a rule)
rule get_rrna:
    output:
        fasta = "rRNA-DB/rrna.fasta",
        idx="rRNA-DB/rrna.fasta.amb"
    conda: "align_env"
    shell:
        """
        (mkdir -p rRNA-DB && \
         cd rRNA-DB && \
         wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz && \
         gunzip -c Homo_sapiens.GRCh38.ncrna.fa.gz | grep "gene_biotype:rRNA" -A 1 --no-group-separator > rrna.fasta && \
         bwa index rrna.fasta && \
         cd ../)
        """


# 4. Trimming raw fastqs
rule fastp:
    input:
        r1="fastq-raw/{sample}_R1.fastq.gz",
        r2="fastq-raw/{sample}_R2.fastq.gz"
    output:
        r1="fastq-trim/{sample}_R1.fastq.gz",
        r2="fastq-trim/{sample}_R2.fastq.gz",
        json="fastq-trim/{sample}.json",
        html="fastq-trim/{sample}.html"
    threads: 8
    conda: "qc_env"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
            -j {output.json} -h {output.html} \
            --detect_adapter_for_pe --trim_poly_g --cut_front --cut_tail \
            --cut_window_size 4 --cut_mean_quality 20 \
            --length_required 36 -w {threads}
        """

# 5. Filter out rRNA (Decontamination)
rule filter_rrna:
    input:
        r1="fastq-trim/{sample}_R1.fastq.gz",
        r2="fastq-trim/{sample}_R2.fastq.gz",
        ref="rRNA-DB/rrna.fasta",
        idx="rRNA-DB/rrna.fasta.amb"
    output:
        r1="fastq-final/{sample}_R1.fastq.gz",
        r2="fastq-final/{sample}_R2.fastq.gz"
    threads: 20
    conda: "align_env"
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | \
        samtools fastq -f 12 -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n
        """

# 6. FastQC on final cleaned files
rule fastqc:
    input:
        "fastq-final/{sample}_{read}.fastq.gz"
    output:
        "fastqc_2/{sample}_{read}_fastqc.html"
    conda: "qc_env"
    shell:
        "fastqc {input} -o fastqc_2/"

# 7. MultiQC (Combines all FastQC reports)
rule multiqc:
    input:
        expand("fastqc_2/{sample}_{read}_fastqc.html", sample=SAMPLES, read=['R1', 'R2'])
    output:
        "fastqc_2.html"
    conda: "qc_env"
    shell:
        "multiqc fastqc_2/ -o ./ -n fastqc_2"

# 8. Download reference transcriptome
rule get_transcriptome:
    output:
        "ref_transcriptome/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    shell:
        """
        (mkdir -p ref_transcriptome && \
         cd ref_transcriptome && \
         wget -c ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz && \
         cd ../)
        """


# 9. Index the transcriptome for Salmon
rule index:
    input:
        ref = "ref_transcriptome/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    output:
        idx = directory("ref_transcriptome/GRCh38_idx")
    threads: 20
    conda: "align_env"
    shell:
        "salmon index -t {input.ref} -i {output.idx} --threads {threads} --keepDuplicates"

# 10. Quantify transcripts
rule quasialign:
    input:
        idx="ref_transcriptome/GRCh38_idx",
        r1 = "fastq-final/{sample}_R1.fastq.gz",
        r2 = "fastq-final/{sample}_R2.fastq.gz"
    output:
        out = "quant/{sample}/quant.sf"
    params:
        dir = "quant/{sample}"
    threads: 10
    conda: "align_env"
    shell:
        "salmon quant -l A -i {input.idx} -1 {input.r1} -2 {input.r2} -o {params.dir} -p {threads} --validateMappings"

