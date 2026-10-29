#!/bin/bash
# =============================================================================
# 06a_blast_asvs.sh -- BLAST every bacterial ASV against NCBI 16S RefSeq (type
# strains). The paper named Erwinia amylovora ASVs and the DMM drivers via BLAST.
#
# INPUT : results/02_qiime2/export/dna-sequences.fasta, reference/blast16S/
# OUTPUT: results/03_R/blast_16S_hits.tsv  (top 5 hits per ASV, outfmt 6 + title)
# =============================================================================
set -euo pipefail
CC=/projects/intro2gds/Team3/Noah/Claude_Code
BIN=$CC/software/envs/sra_tools/bin
mkdir -p "$CC/results/03_R"
"$BIN/blastn" -task megablast \
    -query "$CC/results/02_qiime2/export/dna-sequences.fasta" \
    -db "$CC/reference/blast16S/16S_ribosomal_RNA" \
    -max_target_seqs 5 -num_threads "${SLURM_CPUS_PER_TASK:-4}" \
    -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore stitle" \
    > "$CC/results/03_R/blast_16S_hits.tsv"
echo "[blast] $(cut -f1 "$CC/results/03_R/blast_16S_hits.tsv" | sort -u | wc -l) ASVs with hits"
