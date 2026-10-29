#!/bin/bash
# =============================================================================
# 02c_restart_test.sh -- "What if the download job stopped after 30 samples?"
# Simulates that failure in a separate folder (never touches raw_data/fastq):
#   * copies only the first 30 finished runs (as if the job died there)
#   * leaves a half-written *.tmp file for run 31 (killed mid-write)
#   * truncates one of the 30 files (silent corruption, e.g. disk/quota problem)
# then (1) detects the problem with 02b_check_downloads.sh, (2) resubmits the
# SAME Slurm script, and (3) re-checks + compares md5 with the original download.
#
# INPUT : raw_data/fastq/ (completed download)
# OUTPUT: raw_data/restart_test/restart_test_report.txt  (+ download_report.tsv)
# Usage : bash scripts/02c_restart_test.sh      (waits for the Slurm job)
# =============================================================================
set -uo pipefail
CC=/projects/intro2gds/Team3/Noah/Claude_Code
T=$CC/raw_data/restart_test
export FASTQ_DIR=$T/fastq
rm -rf "$T"; mkdir -p "$FASTQ_DIR"
LOG=$T/restart_test_report.txt
exec > >(tee "$LOG") 2>&1

mapfile -t RUNS < <(python3 -c "
import csv
for r in csv.DictReader(open('$CC/metadata/SraRunInfo.csv')): print(r['Run'])")
echo "== simulate: job killed after ${#RUNS[@]} -> 30 runs finished"
for r in "${RUNS[@]:0:30}"; do cp "$CC/raw_data/fastq/$r.fastq.gz" "$CC/raw_data/fastq/$r.done" "$FASTQ_DIR/"; done
head -c 1000000 "$CC/raw_data/fastq/${RUNS[30]}.fastq.gz" > "$FASTQ_DIR/${RUNS[30]}.fastq.gz.tmp"   # half-written
truncate -s 500000 "$FASTQ_DIR/${RUNS[5]}.fastq.gz"                                                 # corrupted
echo "   copied 30 runs, left ${RUNS[30]}.fastq.gz.tmp (partial), truncated ${RUNS[5]}.fastq.gz"

echo; echo "== 1. detect"
bash "$CC/scripts/02b_check_downloads.sh"; echo "   check exit code: $? (non-zero = incomplete)"
awk -F'\t' 'NR>1 && $5!="OK"{print "   "$1, $5}' "$T/download_report.tsv" | sort | uniq -c -f1 | head -3
echo "   MISSING: $(grep -c MISSING "$T/download_report.tsv")  BAD: $(grep -cP '\tBAD$' "$T/download_report.tsv")"

echo; echo "== 2. restart: resubmit the same script (FASTQ_DIR exported)"
JID=$(sbatch --parsable --export=ALL,FASTQ_DIR="$FASTQ_DIR" --job-name=sra_restart_test \
      --output="$CC/logs/02c_restart_test_%j.out" "$CC/scripts/02_download_sra.slurm")
echo "   submitted job $JID"
while squeue -h -j "$JID" 2>/dev/null | grep -q .; do sleep 15; done
sacct -j "$JID" -X -n -o JobID,State,Elapsed
echo "   SKIP (already verified): $(grep -c '^SKIP' "$CC/logs/02c_restart_test_$JID.out")"
echo "   OK   (downloaded now)  : $(grep -c '^OK' "$CC/logs/02c_restart_test_$JID.out")"
grep -E "^OK +${RUNS[5]}" "$CC/logs/02c_restart_test_$JID.out" | sed 's/^/   corrupted run re-fetched: /'

echo; echo "== 3. verify"
bash "$CC/scripts/02b_check_downloads.sh"; echo "   check exit code: $?"
ls "$FASTQ_DIR"/*.tmp 2>/dev/null | sed 's/^/   leftover tmp: /' || true
# .fastq.gz bytes differ between downloads (gzip stores a timestamp), so compare
# md5 of the decompressed FASTQ content with the original download
nd=0
for r in "${RUNS[@]}"; do
    a=$(zcat "$FASTQ_DIR/$r.fastq.gz" | md5sum | cut -d' ' -f1)
    b=$(zcat "$CC/raw_data/fastq/$r.fastq.gz" | md5sum | cut -d' ' -f1)
    [[ $a == "$b" ]] || { nd=$((nd+1)); echo "   content differs: $r"; }
done
echo "   FASTQ content identical to original download for $(( ${#RUNS[@]} - nd ))/${#RUNS[@]} runs"
rm -f "$FASTQ_DIR"/*.fastq.gz          # free 2 GB; keep reports + .done markers
