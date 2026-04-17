#!/bin/bash
#SBATCH --job-name=setup
#SBATCH --output=icca/setup.out
#SBATCH --error=icca/setup.err
#SBATCH --nodes=1
#SBATCH --mem=20G


# 1. Create the directory first
mkdir -p ~/icca
cd icca

# 2. Load Anaconda
module load anaconda/2024.10

# 3. Create the environments
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n align_env bwa samtools salmon
conda create -n qc_env fastqc multiqc fastp
conda create -n snakemake_env python=3.10 snakemake fastq-dl


