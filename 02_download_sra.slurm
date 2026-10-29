#!/bin/bash
# =============================================================================
# 02_download_sra.slurm  --  download all 48 PacBio HiFi runs from NCBI SRA
#
# INPUT : metadata/SraRunInfo.csv   (list of SRR accessions + expected read counts)
# OUTPUT: raw_data/fastq/<SRR>.fastq.gz   one gzipped FASTQ per run
#         raw_data/fastq/<SRR>.done       marker = run finished + verified
#         raw_data/MD5SUMS                md5 of every .fastq.gz (re-check: md5sum -c)
#         raw_data/download_report.tsv    per-run: expected vs observed reads, status
#
# Verification per run:
#   1. prefetch   - downloads the .sra file; checks it against NCBI's md5
#   2. vdb-validate - re-checks the md5/CRC checksums stored inside the .sra file
#   3. fasterq-dump - converts .sra -> FASTQ; read count must equal RunInfo "spots"
#   4. md5sum of the final .fastq.gz is written to raw_data/MD5SUMS
#
# Restart-safe: a run is skipped if <SRR>.done exists AND its md5 still matches.
# If the job dies after 30 runs, just `sbatch` it again - only the missing /
# corrupted runs are re-downloaded (partial files are written as *.tmp first).
# Check status any time with:  bash scripts/02b_check_downloads.sh
# =============================================================================
#SBATCH --job-name=sra_download          # name shown in squeue
#SBATCH --account=introtogds             # allocation charged for the job
#SBATCH --partition=normal_q             # CPU partition on Tinkercliffs
#SBATCH --nodes=1                        # one node
#SBATCH --ntasks=1                       # one task (the script)...
#SBATCH --cpus-per-task=8                # ...with 8 cores -> 8 parallel downloads
#SBATCH --mem=32G                        # RAM for 8 concurrent fasterq-dump
#SBATCH --time=02:00:00                  # wall-time limit (hh:mm:ss)
#SBATCH --output=/projects/intro2gds/Team3/Noah/Claude_Code/logs/02_download_%j.out   # stdout+stderr (%j = job id)

set -euo pipefail
CC=/projects/intro2gds/Team3/Noah/Claude_Code
BIN=$CC/software/envs/sra_tools/bin           # sra-tools installed by 00_install_software.sh
OUT=${FASTQ_DIR:-$CC/raw_data/fastq}     # override with FASTQ_DIR=... (used by the restart test)
SRA=$CC/raw_data/sra_cache
TMP=${TMPDIR:-/tmp}/sra_$SLURM_JOB_ID
mkdir -p "$OUT" "$SRA" "$TMP"
export PATH=$BIN:$PATH
NJ=${SLURM_CPUS_PER_TASK:-4}

echo "job $SLURM_JOB_ID on $(hostname)  start $(date)"

# list of "SRR expected_spots" from the run table
python3 -c "
import csv
for r in csv.DictReader(open('$CC/metadata/SraRunInfo.csv')): print(r['Run'], r['spots'])
" > "$TMP/runs.txt"

download_one() {
    set -euo pipefail
    local run=$1 spots=$2
    local fq=$OUT/$run.fastq.gz
    # ---- skip if already done and still intact --------------------------------
    if [[ -f $OUT/$run.done && -f $fq ]] && (cd "$OUT" && md5sum -c --status "$run.done"); then
        echo "SKIP  $run (already verified)"; return 0
    fi
    rm -f "$fq" "$fq.tmp" "$OUT/$run.done"
    # ---- 1. download .sra (prefetch verifies NCBI md5) ------------------------
    prefetch --max-size 50G --output-directory "$SRA" "$run" >  "$TMP/$run.log" 2>&1
    # ---- 2. integrity check of the .sra container ------------------------------
    vdb-validate "$SRA/$run" >> "$TMP/$run.log" 2>&1
    # ---- 3. convert to FASTQ ---------------------------------------------------
    fasterq-dump --threads 2 --temp "$TMP" --outdir "$TMP/$run" "$SRA/$run" >> "$TMP/$run.log" 2>&1
    local f=$(ls "$TMP/$run"/*.fastq)
    local n=$(( $(wc -l < "$f") / 4 ))
    if [[ $n -ne $spots ]]; then
        echo "FAIL  $run reads=$n expected=$spots"; return 1
    fi
    pigz -p 2 -c "$f" > "$fq.tmp" && mv "$fq.tmp" "$fq"      # atomic: no half-written .fastq.gz
    rm -rf "$TMP/$run" "$SRA/$run"                           # free space (.sra no longer needed)
    # ---- 4. record md5 of final file; the .done marker IS the md5 line ---------
    (cd "$OUT" && md5sum "$run.fastq.gz" > "$run.done")
    echo "OK    $run reads=$n"
}
export -f download_one
export OUT SRA TMP

# run NJ downloads in parallel; collect failures instead of aborting the batch
set +e
xargs -P "$NJ" -L 1 bash -c 'download_one "$0" "$1"' < "$TMP/runs.txt" | tee "$TMP/status.txt"
set -e

cat "$OUT"/*.done | sort -k2 > "$(dirname "$OUT")/MD5SUMS"
bash "$CC/scripts/02b_check_downloads.sh"
echo "end $(date)"
