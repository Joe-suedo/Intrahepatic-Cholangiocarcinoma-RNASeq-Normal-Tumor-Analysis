#!/bin/bash
#SBATCH --job-name=qc_quasialign
#SBATCH --output=qc_quasialign.out
#SBATCH --error=qc_quasialign.err
#SBATCH --nodes=1
#SBATCH --mem=90G
#SBATCH --cpus-per-task=40

#1 Load Anaconda
module load anaconda/2024.10

#2 Initialize bash
source $(conda info --base)/etc/profile.d/conda.sh

#3 Clear PYTHONPATH to avoid library conflicts
export PYTHONPATH=""

#4 Activate snakemake env
conda activate snakemake_env

#5 Run snakemake file
snakemake --use-conda --cores 20 -p -s qc_quasialign.py
# -p: prints the shell commands (helpful for debugging)
# -s: points to your specific Snakefile
