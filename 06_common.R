# 06_common.R -- shared loader, sourced by 06_community_analysis.R and 07_network_analysis.R
# Reads the QIIME2 exports + metadata + BLAST hits and builds:
#   ps      phyloseq object, raw counts (bacteria only), ASVs renamed ASV0001.. by abundance
#   ps_rar  rarefied copy (paper: "ASV abundances were rarefied to the same sequencing depth")
#   ea_asvs character vector of ASVs identified as Erwinia amylovora
suppressPackageStartupMessages({
  library(phyloseq); library(vegan); library(ggplot2); library(dplyr); library(tidyr)
  library(reshape2); library(magrittr); library(ape); library(Biostrings)
})
CC  <- Sys.getenv("CC", "/projects/intro2gds/Team3/Noah/Claude_Code")
EX  <- file.path(CC, "results/02_qiime2/export")
RD  <- file.path(CC, "results/03_R"); dir.create(RD, FALSE, TRUE)
PDF <- function(name) file.path(CC, name)                     # all PDFs go to Claude_Code/
set.seed(2026)

CULT_ORDER <- c("Cripps Pink", "Aztec Fuji", "Gala", "Hewe's Crab", "Granniwinkle", "Unknown cultivar")
CULT_COL   <- c("Aztec Fuji" = "#F8766D", "Cripps Pink" = "#B79F00", "Gala" = "#00BA38",
                "Granniwinkle" = "#00BFC4", "Hewe's Crab" = "#619CFF", "Unknown cultivar" = "#F564E3")
HS_COL     <- c(Canker = "#c0392b", Healthy = "#27ae60")

# ---- metadata ----------------------------------------------------------------
md <- read.delim(file.path(CC, "metadata/sample_metadata.tsv"), check.names = FALSE, comment.char = "")
md <- md[md[[1]] != "#q2:types", ]
rownames(md) <- md[["sample-id"]]
md$Cultivar      <- factor(md$Cultivar, CULT_ORDER)
md$Health_status <- factor(md$Health_status, c("Canker", "Healthy"))

# ---- feature table / taxonomy / tree ----------------------------------------
read_biom_tsv <- function(f) {
  x <- read.delim(f, skip = 1, check.names = FALSE, row.names = 1)
  as.matrix(x)
}
otu  <- read_biom_tsv(file.path(EX, "feature-table.tsv"))
taxr <- read.delim(file.path(EX, "taxonomy.tsv"), check.names = FALSE)
RANKS <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
split_tax <- function(s) {
  v <- trimws(strsplit(s, ";")[[1]]); v <- sub("^[a-z]__", "", v)
  v <- c(v, rep("", 7 - length(v)))[1:7]; v[v == ""] <- NA; v
}
tax <- t(sapply(taxr$Taxon, split_tax)); dimnames(tax) <- list(taxr[["Feature ID"]], RANKS)
tax <- tax[rownames(otu), , drop = FALSE]
tree <- read.tree(file.path(EX, "tree.nwk"))
seqs <- readDNAStringSet(file.path(EX, "dna-sequences.fasta"))

# rename hashed feature IDs -> ASV0001 ... (ordered by total abundance), keep the map
ord <- order(rowSums(otu), decreasing = TRUE)
newid <- setNames(sprintf("ASV%04d", seq_along(ord)), rownames(otu)[ord])
write.table(data.frame(ASV = newid, Feature_ID = names(newid), Total_reads = rowSums(otu)[names(newid)],
                       tax[names(newid), ]), file.path(RD, "asv_id_map.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
rownames(otu) <- newid[rownames(otu)]; rownames(tax) <- newid[rownames(tax)]
tree$tip.label <- newid[tree$tip.label]; names(seqs) <- newid[names(seqs)]

# fill unassigned ranks with the last known rank (e.g. "Comamonadaceae_g__")
tax_f <- tax
for (i in 2:7) { na <- is.na(tax_f[, i]); tax_f[na, i] <- ifelse(grepl("_$", tax_f[na, i - 1]),
                 tax_f[na, i - 1], paste0(tax_f[na, i - 1], "_", tolower(substr(RANKS[i], 1, 1)), "__")) }

ps <- phyloseq(otu_table(otu[, rownames(md)], taxa_are_rows = TRUE), tax_table(tax_f),
               sample_data(md), phy_tree(tree), refseq(seqs))

# ---- BLAST-based species names (paper: species via NCBI BLAST) --------------
bl <- file.path(RD, "blast_16S_hits.tsv")
blast_top <- NULL
if (file.exists(bl) && file.info(bl)$size > 0) {
  b <- read.delim(bl, header = FALSE, quote = "",
                  col.names = c("q", "s", "pident", "len", "qlen", "slen", "evalue", "bits", "title"))
  b$q <- newid[b$q]
  b$species <- sub("^(\\S+ \\S+( pv\\. \\S+| subsp\\. \\S+)?).*", "\\1", b$title)
  b <- b[order(b$q, -b$bits), ]
  blast_top <- b[!duplicated(b$q), c("q", "species", "pident", "len", "qlen", "title")]
  names(blast_top)[1] <- "ASV"
}
# E. amylovora ASVs: SILVA genus Erwinia AND best BLAST hit E. amylovora at >=99% identity
erw <- taxa_names(ps)[tax[, "Genus"] %in% "Erwinia"]
ea_asvs <- if (!is.null(blast_top)) {
  intersect(erw, blast_top$ASV[grepl("Erwinia amylovora", blast_top$species) & blast_top$pident >= 99])
} else erw
label_asv <- function(a) {
  sp <- if (!is.null(blast_top)) blast_top$species[match(a, blast_top$ASV)] else NA
  ifelse(is.na(sp), tax_f[a, "Genus"], sp)
}

# ---- rarefaction --------------------------------------------------------------
depths <- sample_sums(ps)
RARE_DEPTH <- as.numeric(Sys.getenv("RARE_DEPTH", NA))
if (is.na(RARE_DEPTH)) RARE_DEPTH <- min(depths)
ps_rar <- rarefy_even_depth(ps, sample.size = RARE_DEPTH, rngseed = 2026, replace = FALSE,
                            trimOTUs = TRUE, verbose = FALSE)
