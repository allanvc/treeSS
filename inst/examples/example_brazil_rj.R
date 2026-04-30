##=============================================================================
## Example 1: Infant Mortality in Rio de Janeiro, Brazil (2016)
##
## Reproduces Section 5.2 of:
##   Cancado et al. (2025). Environmental and Ecological Statistics, 32, 953-978.
##
## Two analyses are shown:
##   1. The PAPER-FAITHFUL multi-cluster procedure (Sec. 5.1.1):
##      single scan + filter_clusters / get_cluster_regions(n_clusters = N).
##   2. (Optional) iterative_scan with Holm-Bonferroni correction.
##      This is NOT part of the paper - it is a Zhang/Assuncao/Kulldorff (2010)
##      style conditional procedure, kept here as a useful extension.
##=============================================================================

options(bitmapType = "cairo")     # robust PNG rendering on Linux servers

library(treeSS)
library(ggplot2)
library(geobr)
library(sf)


## ---- 1. Load data ----
data(rj_mortality)
data(rj_tree)

cat("=== Rio de Janeiro Infant Mortality 2016 ===\n")
cat("Rows (long):          ", nrow(rj_mortality), "\n")
cat("Unique municipalities:", length(unique(rj_mortality$region_id)), "\n")
cat("Tree nodes:           ", nrow(rj_tree), "\n")
cat("Total deaths:         ", sum(rj_mortality$cases), "\n")


## ---- 2. Run the tree-spatial scan ----
cat("\nRunning tree-spatial scan (nsim=999, n_cores=4)...\n")
system.time({
  result_rj <- treespatial_scan(
    cases       = rj_mortality$cases,
    population  = rj_mortality$live_births,
    region_id   = rj_mortality$region_id,
    x           = rj_mortality$x,
    y           = rj_mortality$y,
    node_id     = rj_mortality$node_id,
    tree        = rj_tree,
    max_pop_pct = 0.50,
    nsim        = 999, seed = 2016,
    n_cores     = 4L
  )
})
print(result_rj)


## ---- 3. Paper-faithful: distinct top clusters via filter_clusters ----
cat("\n=== Distinct top clusters (paper Sec. 5.1.1) ===\n")
fc <- filter_clusters(result_rj)
print(head(fc[, c("node_id", "n_regions", "cases", "expected", "llr", "pvalue")],
            5))


## ---- 4. (Optional) Iterative scan with Holm-Bonferroni ----
## NOT part of Cancado et al. (2025). p-values corrected for multiple testing.
cat("\n=== Iterative scan (extension, with Holm-Bonferroni) ===\n")
iter_rj <- iterative_scan(
  cases       = rj_mortality$cases,
  population  = rj_mortality$live_births,
  region_id   = rj_mortality$region_id,
  x           = rj_mortality$x,
  y           = rj_mortality$y,
  node_id     = rj_mortality$node_id,
  tree        = rj_tree,
  max_iter    = 5, alpha = 0.05,
  nsim        = 999, seed = 2016,
  max_pop_pct = 0.50, n_cores = 4L
)
print(iter_rj)


## ---- 5. Polygons + cluster membership ----
cat("\nDownloading RJ municipal boundaries...\n")
mun <- read_municipality(code_muni = "RJ", year = 2016)
mun$code6 <- as.integer(substr(mun$code_muni, 1, 6))

region_info <- unique(rj_mortality[, c("region_id", "ibge_code", "name")])

cr1   <- merge(get_cluster_regions(result_rj, n_clusters = 1, overlap = FALSE),
               region_info, by = "region_id")
cr2   <- merge(get_cluster_regions(result_rj, n_clusters = 2, overlap = TRUE),
               region_info, by = "region_id")
cr_it <- merge(get_cluster_regions(iter_rj, overlap = TRUE),
               region_info, by = "region_id")


## ---- 6. Plot 1: Most likely cluster ----
mun1 <- merge(mun, cr1, by.x = "code6", by.y = "ibge_code", all.x = TRUE)
mlc  <- result_rj$most_likely_cluster

p1 <- ggplot(mun1) +
  geom_sf(aes(fill = factor(cluster)), color = "gray40",
          linewidth = 0.15, alpha = 0.75) +
  scale_fill_manual(
    values   = c("1" = "#C44E52"),
    na.value = "gray95",
    labels   = c("1" = paste0(mlc$node_id, " (", length(mlc$region_ids),
                              " municipalities)")),
    name     = "Cluster"
  ) +
  labs(
    title    = paste0("Tree-Spatial Scan: most likely cluster (", mlc$node_id, ")"),
    subtitle = paste0("LR=", round(mlc$llr, 1),
                      ", p=", format.pval(result_rj$pvalue, digits = 3),
                      " | ", mlc$cases, " cases (expected ",
                      round(mlc$expected, 1), ")")
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray40"),
        axis.title    = element_blank(),
        axis.text     = element_text(color = "gray50", size = 8))

ggsave("rj_cluster_mlc.png", p1, width = 9, height = 8, dpi = 300)
cat("Saved: rj_cluster_mlc.png\n")


## ---- 7. Plot 2: Top-2 distinct clusters (paper Sec. 5.1.1) ----
mun2 <- merge(mun, cr2, by.x = "code6", by.y = "ibge_code")
palette2 <- c("1" = "#C44E52", "2" = "#4C72B0")

p2 <- ggplot(mun2) +
  geom_sf(aes(fill = factor(cluster)), color = "gray40",
          linewidth = 0.12, alpha = 0.75) +
  scale_fill_manual(values = palette2, na.value = "gray95",
                    name = "Cluster", na.translate = FALSE) +
  facet_wrap(~ panel, nrow = 1) +
  labs(title    = "Tree-Spatial Scan: Top 2 Distinct Clusters",
       subtitle = "Cancado et al. (2025), Sec. 5.1.1") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray40"),
        strip.text    = element_text(face = "bold", size = 10),
        axis.title    = element_blank(),
        axis.text     = element_text(color = "gray50", size = 7),
        legend.position = "none")

ggsave("rj_clusters_top2.png", p2, width = 14, height = 6, dpi = 300)
cat("Saved: rj_clusters_top2.png\n")


## ---- 8. Plot 3 (optional): Iterative scan clusters ----
n_iter <- iter_rj$n_iter
if (n_iter > 0) {
  mun_it <- merge(mun, cr_it, by.x = "code6", by.y = "ibge_code", all.x = TRUE)
  palette_it <- c("#C44E52", "#4C72B0", "#55A868", "#8172B2", "#CCB974")[
    seq_len(n_iter)
  ]
  names(palette_it) <- as.character(seq_len(n_iter))

  n_sig <- sum(iter_rj$clusters$significant, na.rm = TRUE)
  p_it <- ggplot(mun_it) +
    geom_sf(aes(fill = factor(cluster)), color = "gray40",
            linewidth = 0.12, alpha = 0.75) +
    scale_fill_manual(values = palette_it, na.value = "gray95",
                      name = "Iteration", na.translate = FALSE) +
    facet_wrap(~ panel, nrow = 1) +
    labs(title    = paste0("Iterative Scan (extension): ", n_iter,
                            " iterations, ", n_sig,
                            " significant after Holm-Bonferroni"),
         subtitle = "Each panel = one iteration (cases removed before next). NOT in paper.") +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(color = "gray40"),
          strip.text    = element_text(face = "bold", size = 9),
          axis.title    = element_blank(),
          axis.text     = element_text(color = "gray50", size = 7),
          legend.position = "none")

  ggsave("rj_clusters_iterative.png", p_it,
         width = max(8, 4 * n_iter), height = 6, dpi = 300)
  cat("Saved: rj_clusters_iterative.png (", n_iter, " iterations)\n")
}


## ---- 9. Summary ----
cat("\n=== SUMMARY ===\n")
cat("Most likely cluster:\n")
cat("  Node:", mlc$node_id, "(ICD-10)\n")
cat("  Municipalities:", length(mlc$region_ids), "\n")
cat("  Cases:", mlc$cases, "(expected", round(mlc$expected, 1), ")\n")
cat("  Relative risk:", round(mlc$rr, 2), "\n")
cat("  LR:", round(mlc$llr, 2), "  p-value:", result_rj$pvalue, "\n")
cat("\nDistinct clusters via filter_clusters:", nrow(fc), "\n")
cat("Iterative scan: ", n_iter, " iterations, ",
    sum(iter_rj$clusters$significant, na.rm = TRUE),
    " significant after Holm-Bonferroni\n", sep = "")
cat("\nDone!\n")
