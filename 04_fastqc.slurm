#!/bin/bash
# =============================================================================
# 04_fastqc.slurm -- read quality control (paper: "We assessed the quality of the
# amplified sequence reads using FastQC")
#
# INPUT : raw_data/fastq/*.fastq.gz, metadata/sample_metadata.tsv
# OUTPUT: results/01_fastqc/<SRR>_fastqc.{html,zip}   FastQC report per run
#         results/01_fastqc/read_qc_summary.tsv        reads, length, mean Q per run
#         Claude_Code/Fig_S0_read_QC.pdf               length + quality plots (all runs)
# =============================================================================
#SBATCH --job-name=fastqc
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16               # FastQC: one thread per file
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/projects/intro2gds/Team3/Noah/Claude_Code/logs/04_fastqc_%j.out

set -euo pipefail
CC=/projects/intro2gds/Team3/Noah/Claude_Code
OUT=$CC/results/01_fastqc
mkdir -p "$OUT"
bash "$CC/scripts/02b_check_downloads.sh"          # refuse to run on incomplete data

module load FastQC/0.12.1-Java-11
fastqc --threads "$SLURM_CPUS_PER_TASK" --memory 1024 --outdir "$OUT" "$CC"/raw_data/fastq/*.fastq.gz

module reset >/dev/null 2>&1
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
Rscript "$CC/scripts/04b_plot_read_qc.R" "$CC"
