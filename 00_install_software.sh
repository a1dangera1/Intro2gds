#!/bin/bash
# =============================================================================
# 00_install_software.sh
# Installs the tools that are NOT already provided as ARC (Tinkercliffs) modules.
#
#   Already on ARC (loaded with `module load`, nothing installed):
#     - QIIME2/2024.10-amplicon   -> Apptainer image
#                                    /apps/common/software/QIIME2/2024.10/qiime2-amplicon-2024.10.sif
#                                    (same QIIME2 version the paper used)
#     - FastQC/0.12.1-Java-11
#     - R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2 (phyloseq, vegan, randomForest,
#       DirichletMultinomial, igraph, ggplot2, ...)
#
#   Installed here (conda env, via the ARC Miniforge3 module + mamba):
#     - sra-tools  (prefetch / fasterq-dump / vdb-validate)  -> download from SRA
#     - fastspar   (fast C++ SparCC + bootstrap pseudo p-values) -> co-occurrence network
#     - pigz       (parallel gzip)
#     - blast      (BLAST+, identify Erwinia amylovora / DMM-driver ASVs like the paper)
#   Installed here (R packages from CRAN into a private library):
#     - microeco   (FAPROTAX + LEfSe, the package named in the paper)
#
# INPUT : nothing (needs internet)
# OUTPUT: Claude_Code/software/envs/sra_tools   (conda environment)
#         Claude_Code/software/Rlib             (R library)
#         Claude_Code/logs/00_software_versions.txt
#
# Usage:  bash scripts/00_install_software.sh
# =============================================================================
set -euo pipefail

CC=/projects/intro2gds/Team3/Noah/Claude_Code
ENV=$CC/software/envs/sra_tools
RLIB=$CC/software/Rlib
mkdir -p "$CC/software/envs" "$RLIB" "$CC/logs"

# ---- 1. conda environment with sra-tools, fastspar, pigz ---------------------
module load Miniforge3/25.11.0-1
if [[ ! -x $ENV/bin/fasterq-dump ]]; then
    mamba create -y -p "$ENV" -c conda-forge -c bioconda \
        "sra-tools>=3.1" fastspar pigz blast
fi
[[ -x $ENV/bin/blastn ]] || mamba install -y -p "$ENV" -c conda-forge -c bioconda blast

# sra-tools needs a one-time config; point its cache at our project dir so
# prefetch never fills $HOME/ncbi.
mkdir -p "$CC/raw_data/sra_cache"
"$ENV/bin/vdb-config" --set "/repository/user/main/public/root=$CC/raw_data/sra_cache" >/dev/null
"$ENV/bin/vdb-config" --set /LIBS/GUID="$(uuidgen 2>/dev/null || cat /proc/sys/kernel/random/uuid)" >/dev/null

# ---- 2. R packages (CRAN) into a private library -----------------------------
module reset >/dev/null 2>&1   # (purge would also drop the ARC module paths)
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
export R_LIBS_USER=$RLIB
Rscript -e '
  lib <- Sys.getenv("R_LIBS_USER")
  need <- c("microeco")
  miss <- need[!sapply(need, requireNamespace, quietly = TRUE)]
  if (length(miss)) install.packages(miss, lib = lib, repos = "https://cloud.r-project.org")
'

# ---- 3. record versions ------------------------------------------------------
{
  echo "date: $(date)"
  echo "conda env: $ENV"
  "$ENV/bin/fasterq-dump" --version | grep -v "^$"
  "$ENV/bin/blastn" -version | head -1
  "$ENV/bin/fastspar" --version 2>&1 | grep -i "^Version"
  Rscript -e 'cat("R", as.character(getRversion()), "\n"); for (p in c("microeco","phyloseq","vegan","randomForest","DirichletMultinomial","igraph")) cat(p, as.character(packageVersion(p)), "\n")'
  echo "QIIME2 image: /apps/common/software/QIIME2/2024.10/qiime2-amplicon-2024.10.sif"
} | tee "$CC/logs/00_software_versions.txt"
