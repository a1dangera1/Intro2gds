#!/bin/bash
# =============================================================================
# 06_R_analysis.slurm -- statistics, machine learning, functions and networks
#   06a_blast_asvs.sh          BLAST ASVs vs NCBI 16S RefSeq (species names, E. amylovora)
#   06_community_analysis.R    Figs 2-6, S1, S2 (alpha/beta, PERMANOVA, RF, LEfSe, DMM, FAPROTAX)
#   07_network_analysis.R      Figs 7, 8, S3, Table 1 (SparCC/FastSpar networks)
#
# INPUT : results/02_qiime2/export/*, metadata/sample_metadata.tsv
# OUTPUT: Claude_Code/*.pdf, results/03_R/*.tsv, results/04_network/*
# Env   : RARE_DEPTH (reads/sample for rarefaction; default = smallest sample)
# =============================================================================
#SBATCH --job-name=R_analysis
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16               # FastSpar bootstraps run in parallel
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/projects/intro2gds/Team3/Noah/Claude_Code/logs/06_R_analysis_%j.out

set -euo pipefail
export CC=/projects/intro2gds/Team3/Noah/Claude_Code
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
export R_LIBS_USER=$CC/software/Rlib          # microeco installed by 00_install_software.sh
cd "$CC"
bash   scripts/06a_blast_asvs.sh
Rscript scripts/06_community_analysis.R
Rscript scripts/07_network_analysis.R
ls -1 "$CC"/*.pdf
