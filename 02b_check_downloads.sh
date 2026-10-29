#!/bin/bash
# =============================================================================
# 02b_check_downloads.sh  --  is the raw data complete and intact?
#
# INPUT : metadata/SraRunInfo.csv, raw_data/fastq/*.fastq.gz, raw_data/MD5SUMS
# OUTPUT: raw_data/download_report.tsv  (Run, expected reads, observed reads, md5 ok, status)
#         exit code 0 = all 48 runs present & verified, 1 = something to re-run
#
# Usage: bash scripts/02b_check_downloads.sh
#        If anything is MISSING/BAD, resubmit: sbatch scripts/02_download_sra.slurm
# =============================================================================
set -uo pipefail
CC=/projects/intro2gds/Team3/Noah/Claude_Code
OUT=${FASTQ_DIR:-$CC/raw_data/fastq}     # override with FASTQ_DIR=... (used by the restart test)
REP=$(dirname "$OUT")/download_report.tsv
echo -e "Run\tExpected_reads\tObserved_reads\tMD5_ok\tStatus" > "$REP"
bad=0; tot=0; reads=0
while IFS=, read -r run spots; do
    tot=$((tot+1)); fq=$OUT/$run.fastq.gz
    if [[ ! -f $fq || ! -f $OUT/$run.done ]]; then
        echo -e "$run\t$spots\t0\tNA\tMISSING" >> "$REP"; bad=$((bad+1)); continue
    fi
    md5ok=$(cd "$OUT" && md5sum -c --status "$run.done" && echo yes || echo no)
    n=$(( $(zcat "$fq" | wc -l) / 4 ))
    st=OK; [[ $md5ok == yes && $n -eq $spots ]] || { st=BAD; bad=$((bad+1)); }
    reads=$((reads+n))
    echo -e "$run\t$spots\t$n\t$md5ok\t$st" >> "$REP"
done < <(python3 -c "
import csv
for r in csv.DictReader(open('$CC/metadata/SraRunInfo.csv')): print(r['Run']+','+r['spots'])")
echo "[check] $((tot-bad))/$tot runs OK, $reads reads total, $(du -ch "$OUT"/*.fastq.gz 2>/dev/null | tail -1 | cut -f1) FASTQ.gz"
[[ $bad -eq 0 ]]
