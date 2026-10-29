#!/usr/bin/env Rscript
# =============================================================================
# 07_network_analysis.R -- SparCC co-occurrence networks (paper Figs. 7, 8, S3, Table 1)
#
# Paper: SparCC networks built in MENAP; "majority selection" 5 and 0 (canker) and
# 1 (healthy); |SparCC r| >= 0.2; P < 0.01 (SparCC pseudo-P). MENAP is a web
# service, so here SparCC = FastSpar (C++ SparCC) with 1000 bootstraps for the
# pseudo-P, and all MENAP topology metrics are computed with igraph / sna.
# Interpretation of "majority selection" (paper: "0: keep only ASVs ... detected in
# more than half of the samples"):
#   canker  T5 -> ASVs present in >= 5 of the 32 canker samples
#   canker  T0 -> ASVs present in > 50 % of canker samples (majority)
#   healthy T1 -> ASVs present in >= 1 healthy samples with >= 2 reads overall
#
# INPUT : results/02_qiime2/export/* (raw counts, via 06_common.R)
# OUTPUT: results/04_network/<net>/{counts.tsv,cor.tsv,pvalues.tsv,edges.tsv,nodes.tsv}
#         results/04_network/Table1_network_topology.tsv
#         Fig7A_canker_network_T5.pdf, Fig7B_Ea_interactions_T5.pdf,
#         Fig7C_canker_network_T0.pdf, Fig7D_Ea_interactions_T0.pdf,
#         Fig8A_healthy_network_T1.pdf, Fig8B_healthy_top40_ASV_network.pdf,
#         Fig_S3_ZiPi_module_roles.pdf
# =============================================================================
source(file.path(Sys.getenv("CC", "/projects/intro2gds/Team3/Noah/Claude_Code"), "scripts/06_common.R"))
suppressPackageStartupMessages({ library(igraph); library(ggraph); library(patchwork) })
ND  <- file.path(CC, "results/04_network"); dir.create(ND, FALSE, TRUE)
BIN <- file.path(CC, "software/envs/sra_tools/bin")
NT  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
NBOOT <- as.integer(Sys.getenv("NBOOT", "1000"))
R_MIN <- 0.2; P_MAX <- 0.01

counts <- as(otu_table(ps), "matrix")
sets <- list(
  canker_T5  = list(s = rownames(md)[md$Health_status == "Canker"],  rule = function(m) rowSums(m > 0) >= 5),
  canker_T0  = list(s = rownames(md)[md$Health_status == "Canker"],  rule = function(m) rowSums(m > 0) > ncol(m) / 2),
  healthy_T1 = list(s = rownames(md)[md$Health_status == "Healthy"], rule = function(m) rowSums(m > 0) >= 1 & rowSums(m) >= 2))

run_fastspar <- function(name, m) {
  d <- file.path(ND, name); dir.create(file.path(d, "boot"), FALSE, TRUE); dir.create(file.path(d, "bcor"), FALSE, TRUE)
  f <- file.path(d, "counts.tsv")
  write.table(data.frame(`#OTU ID` = rownames(m), m, check.names = FALSE), f, sep = "\t", quote = FALSE, row.names = FALSE)
  sh <- function(...) { cmd <- paste0(...); if (system(cmd) != 0) stop("failed: ", cmd) }
  if (!file.exists(file.path(d, "pvalues.tsv"))) {
    sh(BIN, "/fastspar --otu_table ", f, " --correlation ", d, "/cor.tsv --covariance ", d, "/cov.tsv ",
       "--iterations 50 --exclude_iterations 10 --threshold 0.1 --threads ", NT, " --seed 1 > ", d, "/fastspar.log 2>&1")
    sh(BIN, "/fastspar_bootstrap --otu_table ", f, " --number ", NBOOT, " --prefix ", d, "/boot/b --threads ", NT, " --seed 1 > /dev/null")
    sh("ls ", d, "/boot/*.tsv | xargs -P ", NT, " -I{} sh -c '", BIN, "/fastspar --otu_table {} ",
       "--correlation ", d, "/bcor/cor_$(basename {}) --covariance /dev/null --iterations 5 --threads 1 --seed 1 > /dev/null 2>&1'")
    sh(BIN, "/fastspar_pvalues --otu_table ", f, " --correlation ", d, "/cor.tsv --prefix ", d, "/bcor/cor_b_ ",
       "--permutations ", NBOOT, " --outfile ", d, "/pvalues.tsv --threads ", NT, " > /dev/null")
    unlink(file.path(d, c("boot", "bcor")), recursive = TRUE)
  }
  rd <- function(x) as.matrix(read.delim(file.path(d, x), row.names = 1, check.names = FALSE))
  list(cor = rd("cor.tsv"), p = rd("pvalues.tsv"), dir = d)
}

# stress centrality: number of shortest paths passing through each node
stress_centrality <- function(g) {
  n <- vcount(g); s <- setNames(numeric(n), V(g)$name)
  for (v in seq_len(n - 1)) {
    sp <- suppressWarnings(all_shortest_paths(g, from = v, to = (v + 1):n, weights = NA)$vpaths)
    for (pth in sp) if (length(pth) > 2) { mid <- pth[-c(1, length(pth))]; s[mid] <- s[mid] + 1 }
  }
  s
}
zipi <- function(g, memb) {
  k <- degree(g); A <- as_adjacency_matrix(g, sparse = FALSE)
  z <- p <- numeric(vcount(g))
  for (i in seq_len(vcount(g))) {
    mi <- memb[i]; kin <- sum(A[i, memb == mi])
    ks <- sapply(V(g)[memb == mi], function(j) sum(A[j, memb == mi]))
    s_k <- sd(ks); z[i] <- if (!is.na(s_k) && s_k > 0) (kin - mean(ks)) / s_k else 0
    kit <- tapply(A[i, ], memb, sum); p[i] <- if (k[i] > 0) 1 - sum((kit / k[i])^2) else 0
  }
  role <- ifelse(z > 2.5 & p > 0.62, "Network hub", ifelse(z > 2.5, "Module hub", ifelse(p > 0.62, "Connector", "Peripheral")))
  data.frame(name = V(g)$name, Zi = z, Pi = p, role = role)
}
topology <- function(g) {
  d <- degree(g); bt <- betweenness(g, weights = NA); st <- stress_centrality(g)
  ev <- eigen_centrality(g, weights = NA)$vector; cl <- suppressWarnings(closeness(g, weights = NA))
  fk <- table(d); pl <- tryCatch(summary(lm(log(as.numeric(fk) / sum(fk)) ~ log(as.numeric(names(fk)))))$r.squared, error = function(e) NA)
  nm <- function(x) { w <- names(which(x == max(x, na.rm = TRUE))); paste(sprintf("%s-%s", w, label_asv(w)), collapse = "; ") }
  sg <- sna::as.sociomatrix.sna(as_adjacency_matrix(g, sparse = FALSE))
  c(Nodes = vcount(g), Links = ecount(g), Positive_links = sum(E(g)$r > 0), Negative_links = sum(E(g)$r < 0),
    `R square of power-law` = pl, `Average degree (avgK)` = mean(d),
    `Average clustering coefficient (avgCC)` = transitivity(g, "average", isolates = "zero"),
    `Average path distance (GD)` = mean_distance(g, weights = NA),
    `Geodesic efficiency (E)` = global_efficiency(g, weights = NA),
    `Harmonic geodesic distance (HD)` = 1 / global_efficiency(g, weights = NA),
    `Maximal degree` = max(d), `Nodes with maximal degree` = nm(d),
    `Centralization of degree (CD)` = centr_degree(g)$centralization,
    `Maximal betweenness` = max(bt), `Nodes with maximal betweenness` = nm(bt),
    `Centralization of betweenness (CB)` = centr_betw(g)$centralization,
    `Maximal stress centrality` = max(st), `Nodes with maximal stress centrality` = nm(st),
    `Centralization of stress centrality (CS)` = sum(max(st) - st) / ((vcount(g) - 1) * max(st)),
    `Maximal eigenvector centrality` = max(ev), `Nodes with maximal eigenvector centrality` = nm(ev),
    `Centralization of eigenvector centrality (CE)` = centr_eigen(g)$centralization,
    `Maximal closeness centrality` = max(cl, na.rm = TRUE), `Nodes with maximal closeness centrality` = nm(cl),
    `Centralization of closeness centrality (CCL)` = suppressWarnings(centr_clo(g)$centralization),
    `Density (D)` = edge_density(g), Reciprocity = reciprocity(g), `Transitivity (Trans)` = transitivity(g, "global"),
    `Connectedness (Con)` = sna::connectedness(sg), Efficiency = sna::efficiency(sg),
    Hierarchy = sna::hierarchy(sg, "krackhardt"), Lubness = sna::lubness(sg),
    Modules = NA, Modularity = NA)
}

genus_of <- function(a) tax_f[a, "Genus"]
nets <- list(); topo <- list(); roles <- list()
for (nm in names(sets)) {
  m <- counts[, sets[[nm]]$s]; m <- m[sets[[nm]]$rule(m), , drop = FALSE]
  m <- m[, colSums(m) > 0, drop = FALSE]                       # SparCC cannot use empty samples
  cat(sprintf("[network] %s: %d ASVs x %d samples\n", nm, nrow(m), ncol(m)))
  fs <- run_fastspar(nm, m)
  ut <- which(upper.tri(fs$cor), arr.ind = TRUE)
  ed <- data.frame(from = rownames(fs$cor)[ut[, 1]], to = colnames(fs$cor)[ut[, 2]],
                   r = fs$cor[ut], p = fs$p[ut])
  ed <- ed[abs(ed$r) >= R_MIN & ed$p < P_MAX, ]
  ed$sign <- ifelse(ed$r > 0, "Positive", "Negative")
  g <- graph_from_data_frame(ed, directed = FALSE)
  g <- delete_vertices(g, degree(g) == 0)
  set.seed(1); cm <- cluster_fast_greedy(simplify(g, edge.attr.comb = "first"), weights = NA)
  V(g)$module <- membership(cm); V(g)$genus <- genus_of(V(g)$name); V(g)$label <- label_asv(V(g)$name)
  V(g)$Ea <- V(g)$name %in% ea_asvs; V(g)$degree <- degree(g)
  t1 <- topology(g); t1["Modules"] <- length(unique(membership(cm))); t1["Modularity"] <- modularity(cm)
  topo[[nm]] <- t1
  zp <- zipi(g, membership(cm)); zp$genus <- genus_of(zp$name); zp$net <- nm; roles[[nm]] <- zp
  write.table(ed, file.path(fs$dir, "edges.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(name = V(g)$name, genus = V(g)$genus, label = V(g)$label, module = V(g)$module,
                         degree = V(g)$degree, Ea = V(g)$Ea, zp[, c("Zi", "Pi", "role")]),
              file.path(fs$dir, "nodes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  nets[[nm]] <- list(g = g, ed = ed, m = m)
}
T1 <- data.frame(Parameter = names(topo[[1]]), sapply(topo, function(x) x[names(topo[[1]])]), check.names = FALSE)
write.table(T1, file.path(ND, "Table1_network_topology.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
R <- bind_rows(roles); write.table(R, file.path(ND, "ZiPi_roles.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(T1, row.names = FALSE)
cat("\nModule hubs / connectors:\n"); print(R[R$role != "Peripheral", ], row.names = FALSE)

# ---- plotting --------------------------------------------------------------------
plot_net <- function(g, title, label_all = TRUE) {
  topg <- names(sort(table(V(g)$genus), decreasing = TRUE))[1:min(11, length(unique(V(g)$genus)))]
  V(g)$col <- ifelse(V(g)$genus %in% topg, V(g)$genus, "Other")
  pal <- setNames(c(RColorBrewer::brewer.pal(11, "Paired")[seq_along(topg)], "grey70"), c(topg, "Other"))
  set.seed(3)
  ggraph(g, layout = "fr") +
    geom_edge_arc(aes(colour = sign, width = abs(r)), strength = .15, alpha = .35) +
    scale_edge_colour_manual(values = c(Positive = "#2a9d8f", Negative = "#e76f51")) +
    scale_edge_width(range = c(.2, 1.2), name = "|SparCC r|") +
    geom_node_point(aes(size = degree, fill = col), shape = 21, colour = ifelse(V(g)$Ea, "black", "grey30"),
                    stroke = ifelse(V(g)$Ea, 1.5, .3)) +
    scale_fill_manual(values = pal, name = "Genus") + scale_size(range = c(2, 9)) +
    { if (label_all) geom_node_text(aes(label = ifelse(Ea, paste0("E. amylovora ", name), genus)),
                                    size = 2.4, repel = TRUE, max.overlaps = 30) } +
    theme_void() + labs(title = title,
      subtitle = sprintf("%d nodes, %d links (%d +, %d -); |r|>=%.1f, pseudo-P<%.2f; black ring = E. amylovora ASV",
                         vcount(g), ecount(g), sum(E(g)$r > 0), sum(E(g)$r < 0), R_MIN, P_MAX))
}
ea_bubble <- function(ed, title) {
  e <- ed[ed$from %in% ea_asvs | ed$to %in% ea_asvs, ]
  if (!nrow(e)) return(ggplot() + labs(title = paste(title, "- no significant E. amylovora links")) + theme_void())
  sw <- !(e$from %in% ea_asvs); tmp <- e$from[sw]; e$from[sw] <- e$to[sw]; e$to[sw] <- tmp
  e$EaLab <- paste0("E. amylovora ", e$from)
  e$partner <- ifelse(e$to %in% ea_asvs, paste0("E. amylovora ", e$to), paste0(label_asv(e$to), " (", e$to, ")"))
  write.table(e, sub(".pdf$", ".tsv", file.path(ND, gsub(" ", "_", title))), sep = "\t", quote = FALSE, row.names = FALSE)
  ggplot(e, aes(partner, EaLab, colour = sign, size = abs(r))) + geom_point(alpha = .8) +
    scale_colour_manual(values = c(Positive = "#2a9d8f", Negative = "#e76f51"), name = "Interaction") +
    scale_size(range = c(2, 9), name = "|Correlation|") + theme_bw() +
    theme(axis.text.x = element_text(angle = 50, hjust = 1, face = "italic", size = 7)) +
    labs(x = NULL, y = NULL, title = title,
         subtitle = sprintf("%d positive, %d negative; r range %.2f to %.2f", sum(e$r > 0), sum(e$r < 0), min(e$r), max(e$r)))
}
ggsave(PDF("Fig7A_canker_network_T5.pdf"), plot_net(nets$canker_T5$g, "Fig. 7A  Canker network (majority selection T=5)"), width = 15, height = 12)
ggsave(PDF("Fig7B_Ea_interactions_T5.pdf"), ea_bubble(nets$canker_T5$ed, "Fig. 7B  E. amylovora ASV interactions, canker T=5"), width = 14, height = 5)
ggsave(PDF("Fig7C_canker_network_T0.pdf"), plot_net(nets$canker_T0$g, "Fig. 7C  Canker network (majority selection T=0)"), width = 11, height = 9)
ggsave(PDF("Fig7D_Ea_interactions_T0.pdf"), ea_bubble(nets$canker_T0$ed, "Fig. 7D  E. amylovora ASV interactions, canker T=0"), width = 11, height = 5)
ggsave(PDF("Fig8A_healthy_network_T1.pdf"), plot_net(nets$healthy_T1$g, "Fig. 8A  Healthy-bark network (majority selection T=1)"), width = 14, height = 12)

# Fig 8B: top-40 ASVs (by degree) of the healthy network on a circle
gh <- nets$healthy_T1$g; top40 <- names(sort(degree(gh), decreasing = TRUE))[1:min(40, vcount(gh))]
g40 <- induced_subgraph(gh, top40)
p8b <- ggraph(g40, layout = "circle") + geom_edge_link(aes(colour = sign, width = abs(r)), alpha = .5) +
  scale_edge_colour_manual(values = c(Positive = "#2a9d8f", Negative = "#e76f51")) + scale_edge_width(range = c(.2, 1.5)) +
  geom_node_point(aes(size = degree), colour = "grey20") +
  geom_node_text(aes(label = paste(name, label)), size = 2.3, repel = TRUE) + theme_void() +
  labs(title = "Fig. 8B  Healthy-bark ASV interaction network (top 40 ASVs by degree)",
       subtitle = sprintf("%d significant interactions among these 40 ASVs", ecount(g40)))
ggsave(PDF("Fig8B_healthy_top40_ASV_network.pdf"), p8b, width = 12, height = 11)

pz <- ggplot(R, aes(Pi, Zi, colour = role)) + geom_point(alpha = .7) +
  geom_hline(yintercept = 2.5, linetype = 2) + geom_vline(xintercept = .62, linetype = 2) +
  ggrepel::geom_text_repel(data = R[R$role != "Peripheral", ], aes(label = paste(name, genus)), size = 2.6, show.legend = FALSE) +
  facet_wrap(~net) + theme_bw() +
  labs(x = "Among-module connectivity (Pi)", y = "Within-module connectivity (Zi)",
       title = "Fig. S3  Zi-Pi node roles (module hub Zi>2.5; connector Pi>0.62; network hub both)")
ggsave(PDF("Fig_S3_ZiPi_module_roles.pdf"), pz, width = 14, height = 5)
