#!/bin/bash
# =============================================================================
# 03_get_reference.sh -- SILVA 138 reference used for taxonomy (paper: "Silva 138
# Reference Taxonomy Database"). Full-length 16S reads -> the full-length
# SILVA 138 99% Naive-Bayes classifier pre-trained for QIIME2 2024.10
# (scikit-learn 1.4.2).
#
# INPUT : internet
# OUTPUT: reference/silva-138-99-nb-classifier.qza
# =============================================================================
set -euo pipefail
CC=/projects/intro2gds/Team3/Noah/Claude_Code
REF=$CC/reference
SIF=/apps/common/software/QIIME2/2024.10/qiime2-amplicon-2024.10.sif
mkdir -p "$REF" "$CC/software/qiime2_home"
F=$REF/silva-138-99-nb-classifier.qza
if [[ ! -s $F ]]; then
    curl -sSfL --retry 5 -o "$F.tmp" \
      https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza
    mv "$F.tmp" "$F"
fi
module load apptainer/1.4.0 >/dev/null 2>&1
apptainer exec --bind "$CC/software/qiime2_home:/home/qiime2" --bind "$CC" "$SIF" \
    qiime tools peek "$F"
md5sum "$F" | tee "$REF/MD5SUMS"

# ---- NCBI 16S RefSeq BLAST database (type-strain 16S; used to name ASVs the way
#      the paper did: "identified through ... Basic Local Alignment Search Tool")
# OUTPUT: reference/blast16S/16S_ribosomal_RNA.*
B=$REF/blast16S; mkdir -p "$B"
if [[ ! -s $B/16S_ribosomal_RNA.nsq && ! -s $B/16S_ribosomal_RNA.00.nsq ]]; then
    curl -sSfL --retry 5 -o "$B/16S.tar.gz"     https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz
    curl -sSfL --retry 5 -o "$B/16S.tar.gz.md5" https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz.md5
    # NCBI publishes an md5 next to the archive -> verify before unpacking
    [[ $(md5sum < "$B/16S.tar.gz" | cut -d' ' -f1) == $(cut -d' ' -f1 "$B/16S.tar.gz.md5") ]] \
        || { echo "16S db md5 mismatch"; exit 1; }
    tar -xzf "$B/16S.tar.gz" -C "$B" && rm "$B/16S.tar.gz"
fi
$CC/software/envs/sra_tools/bin/blastdbcmd -db "$B/16S_ribosomal_RNA" -info | head -4
