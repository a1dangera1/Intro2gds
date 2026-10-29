#!/bin/bash
# =============================================================================
# 05_qiime2_dada2_taxonomy.slurm -- QIIME2 2024.10 processing, as in the paper:
#   "quality control, diversity analysis, and taxonomic classification using the
#    QIIME2 (v2024.10) platform. DADA2 was used within QIIME2 to remove ...
#    sequence noise and to detect and remove chimeras ... Sequences were assigned
#    ... using the Silva 138 Reference Taxonomy Database."
#
# INPUT : raw_data/fastq/*.fastq.gz, metadata/sample_metadata.tsv,
#         reference/silva-138-99-nb-classifier.qza
# OUTPUT: results/02_qiime2/
#           demux.qza/.qzv              imported reads + summary
#           table.qza, rep-seqs.qza     DADA2 ASV table + ASV sequences (all ASVs)
#           denoising-stats.qza         reads kept at each DADA2 step
#           taxonomy.qza                SILVA 138 assignment of every ASV
#           table-bac.qza, rep-seqs-bac.qza  after removing chloroplast /
#                                       mitochondria / eukaryote / unassigned
#           rooted-tree.qza             MAFFT + FastTree (for Faith's PD)
#           taxa-barplot.qzv
#           export/  feature-table.tsv, taxonomy.tsv, dna-sequences.fasta,
#                    tree.nwk, denoising_stats.tsv, *_all.* (pre-filter)
# =============================================================================
#SBATCH --job-name=qiime2_dada2
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16               # DADA2 + MAFFT threads
#SBATCH --mem=80G                        # full-length SILVA NB classifier needs ~40 GB per worker
#SBATCH --time=24:00:00
#SBATCH --output=/projects/intro2gds/Team3/Noah/Claude_Code/logs/05_qiime2_%j.out

set -euo pipefail
CC=/projects/intro2gds/Team3/Noah/Claude_Code
Q=$CC/results/02_qiime2
EX=$Q/export
mkdir -p "$Q" "$EX" "$CC/software/qiime2_home"
NT=${SLURM_CPUS_PER_TASK:-8}
SIF=/apps/common/software/QIIME2/2024.10/qiime2-amplicon-2024.10.sif
export TMPDIR=${TMPDIR:-/tmp}
module load apptainer/1.4.0
# run any qiime command inside the container; bind project + a writable home + tmp
qiime() { apptainer exec --bind "$CC/software/qiime2_home:/home/qiime2" --bind "$CC" --bind "$TMPDIR" \
          --env TMPDIR="$TMPDIR" "$SIF" qiime "$@"; }
step() { echo -e "\n===== $(date +%T)  $*"; }

bash "$CC/scripts/02b_check_downloads.sh"        # refuse to run on incomplete data

# ---- 1. manifest + import --------------------------------------------------------
step "import"
awk -F'\t' -v d="$CC/raw_data/fastq" 'NR==1{print "sample-id\tabsolute-filepath"; next}
     $1!="#q2:types"{print $1"\t"d"/"$2".fastq.gz"}' "$CC/metadata/sample_metadata.tsv" > "$Q/manifest.tsv"
[[ -s $Q/demux.qza ]] || qiime tools import --type 'SampleData[SequencesWithQuality]' \
    --input-format SingleEndFastqManifestPhred33V2 --input-path "$Q/manifest.tsv" --output-path "$Q/demux.qza"
[[ -s $Q/demux.qzv ]] || qiime demux summarize --i-data "$Q/demux.qza" --o-visualization "$Q/demux.qzv"

# ---- 2. DADA2 for PacBio CCS/HiFi --------------------------------------------------
# primers from the paper: 27F AGRGTTYGATYMTGGCTCAG, 1492R RGYTACCTTGTTACGACTT
# (q2-dada2 reverse-complements --p-adapter itself and re-orients reads).
# length window 1000-1600 bp and maxEE 2 follow the DADA2 PacBio full-length 16S workflow.
step "dada2 denoise-ccs"
[[ -s $Q/table.qza ]] || qiime dada2 denoise-ccs \
    --i-demultiplexed-seqs "$Q/demux.qza" \
    --p-front AGRGTTYGATYMTGGCTCAG \
    --p-adapter RGYTACCTTGTTACGACTT \
    --p-min-len 1000 --p-max-len 1600 --p-max-ee 2 \
    --p-chimera-method consensus \
    --p-n-threads "$NT" \
    --o-table "$Q/table.qza" --o-representative-sequences "$Q/rep-seqs.qza" \
    --o-denoising-stats "$Q/denoising-stats.qza" --verbose

# ---- 3. taxonomy: SILVA 138 -------------------------------------------------------
step "classify-sklearn"
[[ -s $Q/taxonomy.qza ]] || qiime feature-classifier classify-sklearn \
    --i-classifier "$CC/reference/silva-138-99-nb-classifier.qza" \
    --i-reads "$Q/rep-seqs.qza" --p-n-jobs 1 --p-reads-per-batch 2000 \
    --o-classification "$Q/taxonomy.qza"

# ---- 4. keep bacteria/archaea only (wood DNA -> many plastid/mito reads) ----------
step "filter host organelles"
qiime taxa filter-table --i-table "$Q/table.qza" --i-taxonomy "$Q/taxonomy.qza" \
    --p-include p__ --p-exclude mitochondria,chloroplast,eukaryota \
    --o-filtered-table "$Q/table-bac.qza"
qiime feature-table filter-seqs --i-data "$Q/rep-seqs.qza" --i-table "$Q/table-bac.qza" \
    --o-filtered-data "$Q/rep-seqs-bac.qza"
qiime taxa barplot --i-table "$Q/table-bac.qza" --i-taxonomy "$Q/taxonomy.qza" \
    --m-metadata-file "$CC/metadata/sample_metadata.tsv" --o-visualization "$Q/taxa-barplot.qzv"

# ---- 5. phylogeny (Faith's PD) -----------------------------------------------------
step "phylogeny"
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences "$Q/rep-seqs-bac.qza" --p-n-threads "$NT" \
    --o-alignment "$Q/aligned.qza" --o-masked-alignment "$Q/masked-aligned.qza" \
    --o-tree "$Q/unrooted-tree.qza" --o-rooted-tree "$Q/rooted-tree.qza"

# ---- 6. export to plain files for R -----------------------------------------------
step "export"
x() { qiime tools export --input-path "$1" --output-path "$EX/tmp" >/dev/null; mv "$EX/tmp/$2" "$EX/$3"; rm -rf "$EX/tmp"; }
x "$Q/table-bac.qza"       feature-table.biom  feature-table.biom
x "$Q/table.qza"           feature-table.biom  feature-table_all.biom
x "$Q/taxonomy.qza"        taxonomy.tsv        taxonomy.tsv
x "$Q/rep-seqs-bac.qza"    dna-sequences.fasta dna-sequences.fasta
x "$Q/rep-seqs.qza"        dna-sequences.fasta dna-sequences_all.fasta
x "$Q/rooted-tree.qza"     tree.nwk            tree.nwk
x "$Q/denoising-stats.qza" stats.tsv           denoising_stats.tsv
for b in feature-table feature-table_all; do
    apptainer exec --bind "$CC" "$SIF" biom convert -i "$EX/$b.biom" -o "$EX/$b.tsv" --to-tsv
done
step "done"
