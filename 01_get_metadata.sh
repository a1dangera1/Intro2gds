#!/bin/bash
# =============================================================================
# 01_get_metadata.sh
# Retrieves the run table for BioProject PRJNA1348523 and builds the sample
# metadata (cultivar, health status, replicate) used by every later step.
#
# INPUT : BioProject accession PRJNA1348523 (internet)
# OUTPUT: metadata/SraRunInfo.csv            NCBI run table (1 row per SRA run)
#         metadata/biosample.xml             BioSample attributes (host, isolation_source)
#         metadata/ena_filereport.tsv        ENA mirror of the same runs (fastq_md5, bytes)
#         metadata/sample_metadata.tsv       QIIME2-format metadata (cleaned)  <- used downstream
#         metadata/metadata_corrections.tsv  every value that was changed, old -> new, why
#
# Usage:  bash scripts/01_get_metadata.sh
# =============================================================================
set -euo pipefail
CC=/projects/intro2gds/Team3/Noah/Claude_Code
BP=PRJNA1348523
OUT=$CC/metadata
mkdir -p "$OUT"

# ---- 1. run table: try EDirect first, fall back to the RunInfo endpoint -----
# esearch sometimes hangs on NCBI's side; `timeout` stops it after 120 s and
# the script then downloads the same table straight from the RunInfo URL.
module load EntrezDirect/24.4.20250710 >/dev/null 2>&1 || true
if timeout 120 bash -c "esearch -db sra -query $BP | efetch -format runinfo" > "$OUT/SraRunInfo.csv.tmp" 2>/dev/null \
   && [[ $(grep -c '^SRR' "$OUT/SraRunInfo.csv.tmp") -gt 0 ]]; then
    echo "[metadata] run table from EDirect (esearch | efetch)"
else
    echo "[metadata] EDirect failed/hung -> using RunInfo endpoint"
    curl -sSf --retry 5 "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=$BP" > "$OUT/SraRunInfo.csv.tmp"
fi
grep -v '^$' "$OUT/SraRunInfo.csv.tmp" > "$OUT/SraRunInfo.csv" && rm -f "$OUT/SraRunInfo.csv.tmp"
echo "[metadata] runs: $(grep -c '^SRR' "$OUT/SraRunInfo.csv")"

# ---- 2. BioSample attributes (host = cultivar, isolation_source = health) ----
IDS=$(python3 -c "import csv;print(','.join(r['BioSample'] for r in csv.DictReader(open('$OUT/SraRunInfo.csv'))))")
curl -sSf --retry 5 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id=$IDS&retmode=xml" > "$OUT/biosample.xml"

# ---- 3. ENA file report (gives an independent FASTQ md5 + read count) --------
curl -sSf --retry 5 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$BP&result=read_run&fields=run_accession,read_count,fastq_ftp,fastq_md5,fastq_bytes&format=tsv" > "$OUT/ena_filereport.tsv"

# ---- 4. build + clean the metadata ------------------------------------------
python3 -I "$CC/scripts/01_build_metadata.py" "$OUT"
