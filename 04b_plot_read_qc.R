#!/usr/bin/env Rscript
# 04b_plot_read_qc.R <Claude_Code dir>
# Parses every FastQC fastqc_data.txt (inside the .zip) and plots read-length
# and per-base quality for all 48 runs, coloured by health status.
# INPUT : results/01_fastqc/*_fastqc.zip, metadata/sample_metadata.tsv
# OUTPUT: results/01_fastqc/read_qc_summary.tsv, Fig_S0_read_QC.pdf
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
cc  <- commandArgs(TRUE)[1]
qcd <- file.path(cc, "results/01_fastqc")
md  <- read.delim(file.path(cc, "metadata/sample_metadata.tsv"), comment.char = "", check.names = FALSE)
md  <- md[md[[1]] != "#q2:types", ]

module_tab <- function(lines, name) {
  s <- grep(paste0("^>>", name), lines); e <- grep("^>>END_MODULE", lines); e <- e[e > s][1]
  x <- lines[(s + 2):(e - 1)]
  do.call(rbind, strsplit(x, "\t"))
}
mid <- function(b) sapply(strsplit(b, "-"), function(v) mean(as.numeric(v)))

len <- qual <- summ <- list()
for (z in list.files(qcd, "_fastqc.zip$", full.names = TRUE)) {
  run <- sub("\\.fastq_fastqc\\.zip$|_fastqc\\.zip$", "", basename(z))
  L <- readLines(unz(z, file.path(sub(".zip$", "", basename(z)), "fastqc_data.txt")))
  bs <- module_tab(L, "Basic Statistics")
  q  <- module_tab(L, "Per base sequence quality")
  l  <- module_tab(L, "Sequence Length Distribution")
  qual[[run]] <- data.frame(Run = run, pos = mid(q[, 1]), mean = as.numeric(q[, 2]))
  len[[run]]  <- data.frame(Run = run, len = mid(l[, 1]), n = as.numeric(l[, 2]))
  sq <- module_tab(L, "Per sequence quality scores")
  summ[[run]] <- data.frame(Run = run,
    Reads = as.numeric(bs[bs[, 1] == "Total Sequences", 2]),
    Length_range = bs[bs[, 1] == "Sequence length", 2],
    Median_length = with(len[[run]], len[which(cumsum(n) >= sum(n) / 2)[1]]),
    Mean_read_Q = round(weighted.mean(as.numeric(sq[, 1]), as.numeric(sq[, 2])), 1))
}
S <- merge(md[, c("sample-id", "Run", "Cultivar", "Health_status")], bind_rows(summ), by = "Run")
write.table(S, file.path(qcd, "read_qc_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

Q <- merge(bind_rows(qual), md[, c("Run", "Health_status")]); Ln <- merge(bind_rows(len), md[, c("Run", "Health_status")])
cols <- c(Canker = "#c0392b", Healthy = "#27ae60")
p1 <- ggplot(Ln[Ln$len > 1000 & Ln$len < 1700, ], aes(len, n, group = Run, colour = Health_status)) +
  geom_line(alpha = .5) + scale_colour_manual(values = cols) + theme_bw() +
  labs(x = "Read length (bp)", y = "Reads", title = "Read-length distribution (PacBio HiFi, 27F-1492R)")
p2 <- ggplot(Q, aes(pos, mean, group = Run, colour = Health_status)) + geom_line(alpha = .5) +
  scale_colour_manual(values = cols) + theme_bw() + coord_cartesian(ylim = c(0, NA)) +
  labs(x = "Position in read (bp)", y = "Mean Phred quality", title = "Per-base quality (FastQC)")
p3 <- ggplot(S, aes(reorder(`sample-id`, Reads), Reads, fill = Health_status)) + geom_col() +
  scale_fill_manual(values = cols) + theme_bw() + coord_flip() +
  labs(x = NULL, y = "Raw HiFi reads", title = sprintf("%d runs, %s reads", nrow(S), format(sum(S$Reads), big.mark = ",")))
ggsave(file.path(cc, "Fig_S0_read_QC.pdf"), (p1 / p2) | p3, width = 15, height = 9)
cat("read QC:", nrow(S), "runs,", sum(S$Reads), "reads; median length", median(S$Median_length), "\n")
