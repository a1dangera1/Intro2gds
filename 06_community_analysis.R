#!/usr/bin/env Rscript
# =============================================================================
# 06_community_analysis.R -- reproduces paper Figs. 2-6 (+ S1, S2)
#
# INPUT : results/02_qiime2/export/*, metadata/sample_metadata.tsv,
#         results/03_R/blast_16S_hits.tsv   (via 06_common.R)
# OUTPUT (PDF, in Claude_Code/):
#   Fig_S1_ASVs_per_sample.pdf, Fig2A_class_relative_abundance.pdf,
#   Fig2B_alpha_diversity.pdf, Fig_S2_richness_canker_vs_healthy.pdf,
#   Fig3_PCoA_BrayCurtis.pdf, Fig4A_heatmap_top20_genera.pdf,
#   Fig4B_RandomForest_cultivar.pdf, Fig5A_LEfSe.pdf, Fig5B_DMM_heatmap.pdf,
#   Fig5C_Erwinia_amylovora_ASVs.pdf, Fig5D_DMM_top10_drivers.pdf,
#   Fig6_FAPROTAX.pdf
# OUTPUT (tables, results/03_R/): alpha_diversity.tsv, alpha_tests.tsv, permanova.tsv,
#   random_forest_importance.tsv, lefse_results.tsv, dmm_fit.tsv, dmm_drivers.tsv,
#   ea_asv_summary.tsv, faprotax_*.tsv, read_tracking.tsv, summary_numbers.tsv
# =============================================================================
source(file.path(Sys.getenv("CC", "/projects/intro2gds/Team3/Noah/Claude_Code"), "scripts/06_common.R"))
suppressPackageStartupMessages({
  library(agricolae); library(randomForest); library(DirichletMultinomial); library(pheatmap)
  library(patchwork); library(microeco); library(RColorBrewer)
})
summ <- list()                                     # headline numbers for the report
add  <- function(k, v) summ[[k]] <<- as.character(v)

# ---- 0. read tracking: raw -> DADA2 -> bacteria only -------------------------
st  <- read.delim(file.path(EX, "denoising_stats.tsv"), check.names = FALSE, comment.char = "#")
all <- read_biom_tsv(file.path(EX, "feature-table_all.tsv"))
trk <- data.frame(sample = st[[1]], raw = st$input, primer_removed = st$`primer-removed`,
                  filtered = st$filtered, denoised = st$denoised, non_chimeric = st$`non-chimeric`)
trk$bacterial <- sample_sums(ps)[trk$sample]
trk <- merge(md[, c("sample-id", "Cultivar", "Health_status")], trk, by.x = "sample-id", by.y = "sample")
write.table(trk, file.path(RD, "read_tracking.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
add("raw_reads", sum(trk$raw)); add("non_chimeric_reads", sum(trk$non_chimeric))
add("ASVs_all_after_chimera_removal", nrow(all)); add("ASVs_bacterial", ntaxa(ps))
add("bacterial_reads", sum(sample_sums(ps))); add("rarefaction_depth", RARE_DEPTH)

# ---- S1: number of ASVs per sample (paper: canker 16-470, mean 130; healthy 1-22, mean 5)
obs <- colSums(otu_table(ps_rar) > 0)
d1 <- data.frame(md[names(obs), ], ASVs = obs)
for (h in c("Canker", "Healthy")) {
  v <- d1$ASVs[d1$Health_status == h]
  add(paste0("ASVs_per_sample_", h), sprintf("range %d-%d, mean %.0f", min(v), max(v), mean(v)))
}
p <- ggplot(d1, aes(reorder(`sample-id`, ASVs), ASVs, fill = Health_status)) + geom_col() +
  facet_grid(. ~ Health_status, scales = "free_x", space = "free") + scale_fill_manual(values = HS_COL) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = .5, size = 7)) +
  labs(x = NULL, y = sprintf("Unique ASVs (rarefied to %d reads)", RARE_DEPTH),
       title = "Fig. S1  Unique ASVs per sample")
ggsave(PDF("Fig_S1_ASVs_per_sample.pdf"), p, width = 12, height = 5)

# ---- Fig 2A: class-level relative abundance -----------------------------------
cls <- tax_glom(transform_sample_counts(ps, function(x) x / sum(x)), "Class", NArm = FALSE)
dc  <- psmelt(cls); dc$Class <- as.character(dc$Class)
topc <- names(sort(tapply(dc$Abundance, dc$Class, sum), decreasing = TRUE))[1:min(21, length(unique(dc$Class)))]
dc$Class[!dc$Class %in% topc] <- "Other"
dc$Class <- factor(dc$Class, rev(c(topc, "Other")))
dc$Sample <- factor(dc$Sample, rownames(md)[order(md$Health_status, md$Cultivar, as.numeric(md$Replicate))])
pal <- colorRampPalette(brewer.pal(12, "Paired"))(length(levels(dc$Class)))
p <- ggplot(dc, aes(Sample, Abundance * 100, fill = Class)) + geom_col(width = .9) +
  facet_grid(. ~ Health_status + Cultivar, scales = "free_x", space = "free",
             labeller = labeller(Cultivar = function(x) gsub(" ", "\n", x))) +
  scale_fill_manual(values = pal) + theme_bw() + scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size = 6), strip.text = element_text(size = 7),
        legend.text = element_text(size = 7)) +
  labs(x = NULL, y = "Relative abundance (%) (Class level)", title = "Fig. 2A  Class-level composition")
ggsave(PDF("Fig2A_class_relative_abundance.pdf"), p, width = 16, height = 7)
cm <- aggregate(Abundance ~ Class + Health_status, dc, mean)
write.table(cm[order(cm$Health_status, -cm$Abundance), ], file.path(RD, "class_mean_relabund.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
gam <- dc$Abundance[dc$Class == "Gammaproteobacteria" & dc$Health_status == "Canker"]
add("max_Gammaproteobacteria_in_canker_pct", round(100 * max(gam), 1))

# ---- Fig 2B + S2: alpha diversity ---------------------------------------------
er <- estimate_richness(ps_rar, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson"))
rownames(er) <- sample_names(ps_rar)
cnt <- as(otu_table(ps_rar), "matrix")
er$Jackknife <- colSums(cnt > 0) + colSums(cnt == 1)              # first-order jackknife (S + f1)
er$Faith_PD  <- picante::pd(t(cnt), phy_tree(ps_rar), include.root = FALSE)$PD
ad <- cbind(md[rownames(er), c("sample-id", "Cultivar", "Health_status")], er)
ad$Faith_PD[is.na(ad$Faith_PD)] <- 0
write.table(ad, file.path(RD, "alpha_diversity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

tests <- list(); letters_l <- list()
for (m in c("Observed", "Chao1", "ACE", "Jackknife", "Shannon", "Simpson", "Faith_PD")) {
  a1 <- anova(lm(ad[[m]] ~ Health_status, ad))
  w  <- wilcox.test(ad[[m]] ~ ad$Health_status, exact = FALSE)
  grp <- interaction(ad$Cultivar, ad$Health_status, sep = "|", drop = TRUE)
  fit <- aov(ad[[m]] ~ grp)
  hsd <- HSD.test(fit, "grp", group = TRUE)$groups
  letters_l[[m]] <- data.frame(metric = m, grp = rownames(hsd), letter = hsd$groups)
  a_c <- anova(lm(ad[[m]] ~ Cultivar, ad[ad$Health_status == "Canker", ]))
  a_h <- anova(lm(ad[[m]] ~ Cultivar, ad[ad$Health_status == "Healthy", ]))
  tests[[m]] <- data.frame(metric = m,
    mean_canker = mean(ad[[m]][ad$Health_status == "Canker"]), mean_healthy = mean(ad[[m]][ad$Health_status == "Healthy"]),
    ANOVA_F_health = a1$`F value`[1], ANOVA_P_health = a1$`Pr(>F)`[1], Wilcoxon_P_health = w$p.value,
    ANOVA_P_cultivar_within_canker = a_c$`Pr(>F)`[1], ANOVA_P_cultivar_within_healthy = a_h$`Pr(>F)`[1])
}
tests <- bind_rows(tests); write.table(tests, file.path(RD, "alpha_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
lets <- bind_rows(letters_l); write.table(lets, file.path(RD, "alpha_tukey_letters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

ad$grp <- paste(ad$Cultivar, ad$Health_status, sep = "|")
alpha_panel <- function(m, ttl) {
  d <- ad; d$val <- d[[m]]
  o <- d %>% group_by(grp, Health_status) %>% summarise(med = median(val), mx = max(val), .groups = "drop") %>%
    arrange(Health_status, desc(med))
  d$grp <- factor(d$grp, o$grp)
  L <- merge(o, lets[lets$metric == m, ], by = "grp")
  ggplot(d, aes(grp, val, colour = Cultivar)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = .15, size = 1) + scale_colour_manual(values = CULT_COL) +
    geom_text(data = L, aes(grp, mx, label = letter), inherit.aes = FALSE, vjust = -.6, size = 3) +
    facet_grid(. ~ Health_status, scales = "free_x", space = "free") +
    scale_x_discrete(labels = function(x) sub("\\|.*", "", x)) +
    scale_y_continuous(expand = expansion(mult = c(.05, .15))) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(x = NULL, y = ttl)
}
p <- (alpha_panel("Simpson", "Simpson") | alpha_panel("Shannon", "Shannon")) /
     (alpha_panel("Observed", "Observed") | alpha_panel("Chao1", "Chao1")) +
  plot_annotation(title = "Fig. 2B  Alpha diversity per cultivar (letters: Tukey HSD across all 12 groups, P<0.05)",
                  subtitle = sprintf("rarefied to %d reads/sample", RARE_DEPTH))
ggsave(PDF("Fig2B_alpha_diversity.pdf"), p, width = 13, height = 10)

ls2 <- melt(ad[, c("sample-id", "Health_status", "Observed", "Chao1", "ACE", "Jackknife", "Faith_PD", "Shannon", "Simpson")],
            id.vars = c("sample-id", "Health_status"))
pl <- tests %>% transmute(variable = metric, lab = sprintf("ANOVA P=%.2g", ANOVA_P_health))
p <- ggplot(ls2, aes(Health_status, value, fill = Health_status)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .15, size = .8) + facet_wrap(~variable, scales = "free_y", nrow = 1) +
  geom_text(data = pl, aes(x = 1.5, y = Inf, label = lab), inherit.aes = FALSE, vjust = 1.5, size = 3) +
  scale_fill_manual(values = HS_COL) + theme_bw() + theme(legend.position = "none") +
  labs(x = NULL, y = NULL, title = "Fig. S2  Richness / diversity: canker vs healthy")
ggsave(PDF("Fig_S2_richness_canker_vs_healthy.pdf"), p, width = 16, height = 4.5)

# ---- Fig 3: PCoA Bray-Curtis + PERMANOVA ----------------------------------------
rel <- transform_sample_counts(ps_rar, function(x) x / sum(x))
bc  <- phyloseq::distance(rel, "bray")
sd  <- data.frame(sample_data(rel))
pm  <- adonis2(bc ~ Health_status * Cultivar, data = sd, permutations = 999, by = "terms")
pm_hs <- adonis2(bc ~ Health_status, data = sd, permutations = 999)
disp <- anova(betadisper(bc, sd$Health_status))
pm_df <- data.frame(term = rownames(pm), pm, check.names = FALSE)
write.table(pm_df, file.path(RD, "permanova.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# pairwise cultivars within cankers
cn  <- subset_samples(rel, Health_status == "Canker"); bcc <- phyloseq::distance(cn, "bray"); sdc <- data.frame(sample_data(cn))
pw <- combn(levels(droplevels(sdc$Cultivar)), 2, function(pr) {
  k <- sdc$Cultivar %in% pr
  r <- adonis2(as.dist(as.matrix(bcc)[k, k]) ~ Cultivar, data = droplevels(sdc[k, ]), permutations = 999)
  data.frame(A = pr[1], B = pr[2], R2 = r$R2[1], F = r$F[1], P = r$`Pr(>F)`[1])
}, simplify = FALSE) %>% bind_rows() %>% mutate(P_BH = p.adjust(P, "BH"))
write.table(pw, file.path(RD, "permanova_pairwise_cultivar_in_canker.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
add("PERMANOVA_health", sprintf("R2=%.3f F=%.2f P=%.3f", pm$R2[1], pm$F[1], pm$`Pr(>F)`[1]))
add("PERMANOVA_cultivar", sprintf("R2=%.3f F=%.2f P=%.3f", pm$R2[2], pm$F[2], pm$`Pr(>F)`[2]))
add("PERMANOVA_health_x_cultivar", sprintf("R2=%.3f F=%.2f P=%.3f", pm$R2[3], pm$F[3], pm$`Pr(>F)`[3]))
add("betadisper_health_P", sprintf("%.3g", disp$`Pr(>F)`[1]))

pc  <- ordinate(rel, "PCoA", distance = bc)
ev  <- round(100 * pc$values$Relative_eig[1:2], 1)
dd  <- data.frame(pc$vectors[, 1:2], sd)
p <- ggplot(dd, aes(Axis.1, Axis.2)) +
  stat_ellipse(aes(group = Health_status), linetype = 2, level = .95) +
  geom_point(aes(colour = Cultivar, shape = Health_status), size = 3) +
  scale_colour_manual(values = CULT_COL) + scale_shape_manual(values = c(Canker = 16, Healthy = 17)) +
  theme_bw() + labs(x = sprintf("PCoA1 [%.1f%%]", ev[1]), y = sprintf("PCoA2 [%.1f%%]", ev[2]),
    title = "Fig. 3  PCoA, Bray-Curtis",
    subtitle = sprintf("PERMANOVA: health R2=%.3f P=%.3f | cultivar R2=%.3f P=%.3f | health x cultivar R2=%.3f P=%.3f",
                       pm$R2[1], pm$`Pr(>F)`[1], pm$R2[2], pm$`Pr(>F)`[2], pm$R2[3], pm$`Pr(>F)`[3]))
ggsave(PDF("Fig3_PCoA_BrayCurtis.pdf"), p, width = 9, height = 7)

# ---- Fig 4A: heatmap of top-20 genera in canker samples ---------------------------
cg  <- tax_glom(transform_sample_counts(subset_samples(ps, Health_status == "Canker"), function(x) 100 * x / sum(x)),
                "Genus", NArm = FALSE)
gm  <- as(otu_table(cg), "matrix"); rownames(gm) <- as.character(tax_table(cg)[, "Genus"])
gm  <- rowsum(gm, rownames(gm))
top20 <- names(sort(rowMeans(gm), decreasing = TRUE))[1:20]
add("top5_genera_canker", paste(top20[1:5], collapse = ", "))
scd <- md[colnames(gm), ]; o <- order(scd$Cultivar, as.numeric(scd$Replicate))
hm  <- log10(gm[top20, o] + 0.01)
ann <- data.frame(Cultivar = scd$Cultivar[o], row.names = colnames(gm)[o])
gaps <- cumsum(table(droplevels(scd$Cultivar[o])))
pdf(PDF("Fig4A_heatmap_top20_genera.pdf"), width = 12, height = 7)
pheatmap(hm, cluster_rows = FALSE, cluster_cols = FALSE, gaps_col = gaps[-length(gaps)],
         annotation_col = ann, annotation_colors = list(Cultivar = CULT_COL[levels(droplevels(scd$Cultivar))]),
         color = colorRampPalette(c("navy", "lightyellow", "firebrick"))(100), border_color = "grey60",
         breaks = seq(-2, 2, length.out = 101), legend_breaks = -2:2,
         legend_labels = c("0.01", "0.1", "1", "10", "100"), fontsize_col = 7,
         main = "Fig. 4A  Top 20 genera in cankers (% relative abundance, log scale)")
dev.off()

# ---- Fig 4B: Random Forest -- which genera predict the cultivar of a canker -------
keepg <- rownames(gm)[rowSums(gm > 0) >= 2]                    # genera present in >= 2 cankers
X <- t(gm[keepg, ]); colnames(X) <- make.names(keepg, unique = TRUE)
y <- droplevels(scd[rownames(X), "Cultivar"])
rf <- randomForest(X, y, ntree = 5000, importance = TRUE)
imp <- data.frame(Genus = keepg, importance(rf), check.names = FALSE)
imp <- imp[order(-imp$MeanDecreaseGini), ]
write.table(imp, file.path(RD, "random_forest_importance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
add("RF_OOB_error_pct", round(100 * tail(rf$err.rate[, "OOB"], 1), 1))
top8 <- head(imp, 8)
# which cultivar each top genus is most abundant in
mab <- sapply(top8$Genus, function(g) tapply(gm[g, rownames(X)], y, mean))
top8$Cultivar_max <- rownames(mab)[apply(mab, 2, which.max)]
add("RF_top8_genera", paste(sprintf("%s(%s)", top8$Genus, top8$Cultivar_max), collapse = ", "))
d4 <- melt(mab); names(d4) <- c("Cultivar", "Genus", "RelAb")
d4$Genus <- factor(d4$Genus, rev(top8$Genus)); top8$Genus <- factor(top8$Genus, rev(top8$Genus))
pa <- ggplot(top8, aes(MeanDecreaseGini, Genus, fill = factor(Cultivar_max, CULT_ORDER))) + geom_col() +
  scale_fill_manual(values = CULT_COL, name = "Most abundant in") + theme_bw() + labs(y = NULL)
pb <- ggplot(d4, aes(RelAb, Genus, fill = Cultivar)) + geom_col(position = "dodge") +
  scale_fill_manual(values = CULT_COL) + theme_bw() + labs(y = NULL, x = "Mean relative abundance (%)")
ggsave(PDF("Fig4B_RandomForest_cultivar.pdf"),
       (pa | pb) + plot_annotation(title = "Fig. 4B  Random Forest: top 8 genera discriminating canker microbiomes by cultivar",
                                   subtitle = sprintf("5000 trees, OOB error %.1f%%", 100 * tail(rf$err.rate[, "OOB"], 1))),
       width = 14, height = 6)

# ---- Fig 5A: LEfSe canker vs healthy (microeco, as in the paper) ------------------
tt <- as.data.frame(tax_f, stringsAsFactors = FALSE)
mt <- microtable$new(otu_table = as.data.frame(as(otu_table(ps), "matrix")),
                     sample_table = data.frame(sample_data(ps)), tax_table = tt, phylo_tree = phy_tree(ps))
mt$tidy_dataset()
lf <- trans_diff$new(dataset = mt, method = "lefse", group = "Health_status", alpha = 0.05,
                     lefse_subgroup = NULL, taxa_level = "Genus", p_adjust_method = "none")
lres <- lf$res_diff; write.table(lres, file.path(RD, "lefse_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
ld <- lres[lres$LDA >= 3.5, ]; ld$Taxa <- sub(".*\\|g__", "", ld$Taxa)
ld <- ld[order(ld$LDA), ]; ld$Taxa <- factor(ld$Taxa, unique(ld$Taxa))
add("LEfSe_n_genera_LDA>=3.5", nrow(ld))
add("LEfSe_top_canker", paste(head(rev(as.character(ld$Taxa[ld$Group == "Canker"])), 5), collapse = ", "))
add("LEfSe_top_healthy", paste(head(rev(as.character(ld$Taxa[ld$Group == "Healthy"])), 5), collapse = ", "))
p <- ggplot(ld, aes(LDA, Taxa, colour = Group, size = LDA)) + geom_point() +
  geom_vline(xintercept = c(2, 4), linetype = 2, colour = "grey60") +
  scale_colour_manual(values = HS_COL, name = "Enriched in") + theme_bw() +
  labs(x = "LDA score (log10)", y = "Bacterial genera", title = "Fig. 5A  LEfSe: canker vs healthy (genus, LDA >= 3.5, P<0.05)")
ggsave(PDF("Fig5A_LEfSe.pdf"), p, width = 9, height = max(5, .22 * nrow(ld) + 2))

# ---- Fig 5B/D: Dirichlet-multinomial mixtures, canker samples (n = 32) -------------
cc_ps <- subset_samples(ps, Health_status == "Canker"); cc_ps <- prune_taxa(taxa_sums(cc_ps) > 0, cc_ps)
cnt_c <- t(as(otu_table(cc_ps), "matrix"))
cnt_c <- cnt_c[, colSums(cnt_c > 0) >= 3, drop = FALSE]                # ASVs seen in >= 3 cankers
fits <- lapply(1:7, function(k) dmn(cnt_c, k, verbose = FALSE, seed = 2026))
gof  <- data.frame(k = 1:7, Laplace = sapply(fits, laplace), AIC = sapply(fits, AIC), BIC = sapply(fits, BIC))
write.table(gof, file.path(RD, "dmm_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
kbest <- which.min(gof$Laplace); K <- 4                                  # paper presents k = 4
add("DMM_best_k_Laplace", kbest)
fk <- fits[[K]]
cl <- apply(mixture(fk), 1, which.max)
add("DMM_k4_cluster_sizes", paste(sprintf("C%d=%d", 1:K, tabulate(cl, K)), collapse = " "))
write.table(data.frame(sample = names(cl), cluster = cl, md[names(cl), c("Cultivar")]),
            file.path(RD, "dmm_k4_assignment.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
p0 <- fitted(fits[[1]], scale = TRUE); pk <- fitted(fk, scale = TRUE)
diff <- rowSums(abs(pk - as.vector(p0))); o <- order(diff, decreasing = TRUE)
drv <- data.frame(ASV = rownames(pk)[o], driver_score = diff[o], cumulative = cumsum(diff[o]) / sum(diff),
                  pk[o, ], Genus_SILVA = tax_f[rownames(pk)[o], "Genus"], BLAST = label_asv(rownames(pk)[o]))
write.table(drv, file.path(RD, "dmm_drivers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
top10 <- head(drv, 10)
add("DMM_top10_drivers", paste(sprintf("%s(%s)", top10$ASV, top10$BLAST), collapse = "; "))
top10$lab <- factor(top10$ASV, rev(top10$ASV))
p <- ggplot(top10, aes(driver_score, lab)) + geom_col(fill = "grey35") +
  geom_text(aes(x = max(driver_score) * 1.02, label = BLAST), hjust = 0, fontface = "italic", size = 3.2) +
  scale_x_continuous(expand = expansion(mult = c(0, .9))) + theme_bw() +
  labs(x = "Driver score (sum |fitted k=4 - fitted k=1|)", y = "ASV",
       title = sprintf("Fig. 5D  Top 10 DMM community drivers (k = %d; best k by Laplace = %d)", K, kbest))
ggsave(PDF("Fig5D_DMM_top10_drivers.pdf"), p, width = 10, height = 5)

# heatmap: samples (rows, grouped by DMM cluster) x top-30 driver ASVs
hx <- cnt_c[order(cl), head(drv$ASV, 30)]
rel_hx <- 100 * hx / rowSums(cnt_c[order(cl), ])
colnames(rel_hx) <- paste0(colnames(rel_hx), " ", label_asv(colnames(rel_hx)))
ra <- data.frame(DMM_cluster = factor(cl[order(cl)]), Cultivar = md[rownames(hx), "Cultivar"], row.names = rownames(hx))
pdf(PDF("Fig5B_DMM_heatmap.pdf"), width = 13, height = 9)
pheatmap(sqrt(rel_hx), cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = ra,
         annotation_colors = list(Cultivar = CULT_COL[levels(droplevels(ra$Cultivar))]),
         gaps_row = cumsum(tabulate(cl, K))[-K], fontsize_col = 7, angle_col = 90,
         color = colorRampPalette(c("#3b6fb6", "#f7f7f7", "#d73027"))(100),
         main = sprintf("Fig. 5B  DMM clusters (k=%d) - sqrt(%% relative abundance) of top 30 driver ASVs; sizes %s",
                        K, paste(tabulate(cl, K), collapse = "/")))
plot(gof$k, gof$Laplace, type = "b", xlab = "Number of Dirichlet components (k)", ylab = "Laplace approximation",
     main = "DMM model fit (lower = better)")
dev.off()

# ---- Fig 5C: Erwinia amylovora ASVs per canker sample ------------------------------
add("Ea_ASVs", paste(ea_asvs, collapse = ","))
rel_all <- 100 * t(t(as(otu_table(ps), "matrix")) / sample_sums(ps))
ea_top <- names(sort(rowSums(rel_all[ea_asvs, , drop = FALSE]), decreasing = TRUE))
show <- head(ea_top, 6)
de <- melt(rel_all[show, rownames(md)[md$Health_status == "Canker"], drop = FALSE])
names(de) <- c("ASV", "sample", "RelAb")
ord_s <- names(sort(colSums(rel_all[ea_asvs, rownames(md)[md$Health_status == "Canker"], drop = FALSE])))
de$sample <- factor(de$sample, ord_s); de$ASV <- factor(de$ASV, rev(show))
pa <- ggplot(de, aes(sample, RelAb, fill = ASV)) + geom_col() + theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size = 7)) +
  labs(x = NULL, y = "Relative abundance (%)", title = "Erwinia amylovora ASVs per canker sample")
es <- data.frame(md[md$Health_status == "Canker", c("sample-id", "Cultivar")],
  Ea_relab = colSums(rel_all[ea_asvs, md$Health_status == "Canker", drop = FALSE]),
  Ea_ASVs_detected = colSums(as(otu_table(ps_rar), "matrix")[intersect(ea_asvs, taxa_names(ps_rar)),
                                                             rownames(md)[md$Health_status == "Canker"], drop = FALSE] > 0),
  Ea_shannon = apply(as(otu_table(ps), "matrix")[ea_asvs, md$Health_status == "Canker", drop = FALSE], 2, vegan::diversity))
ecs <- es %>% group_by(Cultivar) %>% summarise(n = n(), mean_Ea_relab = mean(Ea_relab),
         mean_Ea_ASVs = mean(Ea_ASVs_detected), total_Ea_reads = NA, mean_Ea_shannon = mean(Ea_shannon), .groups = "drop")
ecs$total_Ea_reads <- sapply(as.character(ecs$Cultivar), function(cv)
  sum(as(otu_table(ps), "matrix")[ea_asvs, rownames(md)[md$Cultivar == cv & md$Health_status == "Canker"]]))
ecs <- ecs[order(-ecs$mean_Ea_relab), ]
write.table(es, file.path(RD, "ea_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ecs, file.path(RD, "ea_asv_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
add("Ea_rank_by_cultivar_mean_relab", paste(sprintf("%s(%.1f%%)", ecs$Cultivar, ecs$mean_Ea_relab), collapse = " > "))
add("Ea_in_healthy_max_relab_pct", round(max(colSums(rel_all[ea_asvs, md$Health_status == "Healthy", drop = FALSE])), 2))
pb <- ggplot(es, aes(factor(Cultivar, ecs$Cultivar), Ea_relab, fill = Cultivar)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .15) + scale_fill_manual(values = CULT_COL) + theme_bw() + theme(legend.position = "none",
  axis.text.x = element_text(angle = 30, hjust = 1)) + labs(x = NULL, y = "E. amylovora (% of bacterial reads)")
pc2 <- ggplot(es, aes(factor(Cultivar, ecs$Cultivar), Ea_shannon, fill = Cultivar)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .15) + scale_fill_manual(values = CULT_COL) + theme_bw() + theme(legend.position = "none",
  axis.text.x = element_text(angle = 30, hjust = 1)) + labs(x = NULL, y = "Shannon diversity of Ea ASVs")
ggsave(PDF("Fig5C_Erwinia_amylovora_ASVs.pdf"), pa / (pb | pc2) +
       plot_annotation(title = sprintf("Fig. 5C  Erwinia amylovora: %d ASVs (SILVA Erwinia + BLAST top hit E. amylovora >=99%%)",
                                       length(ea_asvs))), width = 13, height = 10)

# ---- Fig 6: FAPROTAX functional profiles (microeco) --------------------------------
tf <- trans_func$new(mt)
tf$cal_func(prok_database = "FAPROTAX")
tf$cal_func_FR(abundance_weighted = TRUE)
fp <- tf$res_func_FR                                           # samples x functions (%)
write.table(data.frame(sample = rownames(fp), fp, check.names = FALSE), file.path(RD, "faprotax_per_sample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
fl <- tf$func_group_list
cat_of <- function(f) { for (g in names(fl)) for (s in names(fl[[g]])) if (f %in% fl[[g]][[s]]) return(s); "Other" }
keep_f <- colnames(fp)[colSums(fp) > 0]
fcat <- sapply(keep_f, cat_of)
fw <- bind_rows(lapply(keep_f, function(f) {
  w <- suppressWarnings(wilcox.test(fp[md$Health_status[match(rownames(fp), rownames(md))] == "Canker", f],
                                    fp[md$Health_status[match(rownames(fp), rownames(md))] == "Healthy", f], exact = FALSE))
  data.frame(Function = f, Category = fcat[f],
             mean_canker = mean(fp[rownames(fp) %in% rownames(md)[md$Health_status == "Canker"], f]),
             mean_healthy = mean(fp[rownames(fp) %in% rownames(md)[md$Health_status == "Healthy"], f]), P = w$p.value)
})) %>% mutate(P_BH = p.adjust(P, "BH")) %>% arrange(P)
write.table(fw, file.path(RD, "faprotax_canker_vs_healthy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
add("FAPROTAX_sig_BH<0.05", paste(sprintf("%s(%s)", fw$Function[fw$P_BH < .05],
    ifelse(fw$mean_canker[fw$P_BH < .05] > fw$mean_healthy[fw$P_BH < .05], "canker", "healthy")), collapse = ", "))
so <- rownames(md)[order(md$Health_status, md$Cultivar, as.numeric(md$Replicate))]
fo <- keep_f[order(factor(fcat, unique(c("Energy source", "C-cycle", "N-cycle", "S-cycle", unique(fcat)))))]
fm <- t(fp[so, fo]); rownames(fm) <- sprintf("%s%s", fo, ifelse(fw$P_BH[match(fo, fw$Function)] < .05, " *", ""))
pdf(PDF("Fig6_FAPROTAX.pdf"), width = 16, height = max(8, .2 * length(fo) + 3))
pheatmap(fm, cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col = data.frame(Health = md[so, "Health_status"], Cultivar = md[so, "Cultivar"], row.names = so),
         annotation_row = data.frame(Category = fcat[fo], row.names = rownames(fm)),
         annotation_colors = list(Health = HS_COL, Cultivar = CULT_COL),
         gaps_col = sum(md$Health_status == "Canker"), color = colorRampPalette(c("white", "#f4a6b7", "#8b0a32"))(100),
         border_color = "grey85", fontsize_row = 7, fontsize_col = 6,
         main = "Fig. 6  FAPROTAX functional groups (% of reads, abundance weighted); * = Wilcoxon canker vs healthy, BH P<0.05")
dev.off()

# ---- headline numbers -------------------------------------------------------------
sn <- data.frame(item = names(summ), value = unlist(summ))
write.table(sn, file.path(RD, "summary_numbers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(sn, row.names = FALSE, right = FALSE)
