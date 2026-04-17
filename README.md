# Intrahepatic Cholangiocarcinoma (iCCA) RNA-Seq Analysis

This repository contains a complete bioinformatics pipeline for Differential Expression Gene (DEG) analysis and variant profiling of Intrahepatic Cholangiocarcinoma using matched Tumor-Normal samples.

## 🧬 Project Overview
The goal of this study is to identify transcriptionally altered pathways and variant profiles in iCCA patients. The analysis covers the entire workflow from raw data acquisition to functional enrichment.

* **Dataset:** BioProject PRJNA488803
* **Samples:** 30 samples (15 matched Tumor-Normal pairs)
* **Organism:** *Homo sapiens*

## 📁 Repository Structure
* **setup.sh**: Shell script for organizing the environments and directory structure.
* **`icca/`**: Shell and Python scripts used for upstream processing on the HPC.
    * Fastq file download and management.
    * Quality Control (FastQC & MultiQC).
    * Transcriptome quantification using **Salmon**.
* **`quant/`**: Output files from Salmon (quant.sf) used for downstream analysis.
* **`Rplots/`**: Visualizations generated during the analysis (Volcano plots, Heatmaps, GSEA).
* **`DESeq2.R`**: Main R script for differential expression analysis.
* **`degs.txt`**: Final list of significantly differentially expressed genes.

## 🚀 Workflow

### 1. Upstream (HPC)
The upstream pipeline was executed on the ACE High-Performance Computing cluster:
1. **QC:** Raw reads were checked for quality.
2. **Quantification:** Reads were quantified against the human reference transcriptome using `Salmon`.

### 2. Downstream (R/Bioconductor)
The downstream analysis was performed in RStudio:
* **Importing:** Transcript quantification estimates were imported with `tximeta`.
* **Normalization:** Using `DESeq2` to handle library size differences and shrinkage.
* **DEGs:** Identified using a threshold of $|log2FC| > 1$ and $padj < 0.05$.
* **Enrichment:** Functional profiling via Gene Set Enrichment Analysis (GSEA) and KEGG pathway analysis.

## 🛠 Tools Used
**Programming Languages:**
* **Bash:** Shell scripting for job submission on the HPC and command-line automation.
* **Python:** Scripting for data handling and integration within the Snakemake workflow.
* **R:** Statistical computing, Bioconductor-based differential expression, and visualization.

**Upstream Analysis on ACE HPC:**
* **Snakemake:** Scalable and reproducible workflow management.
* **Fastq-dl:** Automated download of raw FASTQ files from the SRA.
* **FastQC & MultiQC:** Raw read quality assessment and aggregated reporting.
* **Fastp:** Quality filtering, adapter trimming, and preprocessing.
* **BWA-MEM & Samtools:** Sequence alignment and filtering to remove rRNA (Decontamination).
* **Salmon:** Transcript-level quantification.

**Downstream Analysis in RStudio:**
* **tximeta:** Import of Salmon transcript quantifications.
* **DESeq2:** Statistical normalization and Differential Expression Analysis.
* **AnnotationDbi & org.Hs.eg.db:** Gene ID mapping and human genome annotation.
* **clusterProfiler:** GO and KEGG functional enrichment analysis (GSEA/ORA).
* **enrichplot:** Visualization of pathway analysis results.
* **EnhancedVolcano:** High-quality Volcano plots for DEG distribution.
* **pheatmap:** Hierarchical clustering and expression heatmaps.

## 📊 Key Results and Findings
**Differential Gene Expression**
* PCA and Sample Distance Heatmaps show a distinct transcriptional signature separating Tumor (T) and Normal (N) tissues.
* Observed a dominance of gene downregulation in tumor samples compared to upregulation.

**Functional Enrichment (GO & KEGG)**
* Significant downregulation of pathways related to xenobiotic metabolism, bile acid secretion, and fatty acid degradation—typical of biliary tract cancers.
* Enrichment in cell cycle progression, DNA replication, and extracellular matrix organization, signaling active tumor proliferation and remodeling.

**Pathway Interpretation (GSEA)**
* Results show strong suppression of hepatic metabolic pathways and significant activation of oncogenic signaling.

---
**Author:** Nyanzi Joseph
* Human Genomics Volunteer
* Infectious Diseases Institute - Makerere
* African Center of Excellence in Bioinformatics and Data-Intensive Sciences (ACE)
* Project Duration: 15/10/2025 - 28/02/2026
