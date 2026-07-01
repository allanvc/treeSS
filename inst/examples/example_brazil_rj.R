##=============================================================================
## Example 1: Infant Mortality in Rio de Janeiro, Brazil (2016)
##
## Reproduces Section 5.2 of:
##   Cancado et al. (2025). Environmental and Ecological Statistics, 32, 953-978.
##
## Two secondary-cluster analyses are shown:
##   1. The PAPER-FAITHFUL multi-cluster procedure (Cancado et al. 2025,
##      Sec. 5.1.1): single scan + filter_clusters /
##      get_cluster_regions(n_clusters = N).
##   2. sequential_scan() -- Zhang, Assuncao & Kulldorff (2010).
##=============================================================================

options(bitmapType = "cairo")     # robust PNG rendering on Linux servers

library(treeSS)
library(ggplot2)
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
    rj_mortality,
    cases       = cases,
    population  = live_births,
    region_id   = region_id,
    x           = x,
    y           = y,
    node_id     = node_id,
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


## ---- 4. Sequential scan (Zhang, Assuncao & Kulldorff, 2010) ----
cat("\n=== Sequential scan (Zhang et al. 2010) ===\n")
seq_rj <- sequential_scan(
  rj_mortality,
  cases       = cases,
  population  = live_births,
  region_id   = region_id,
  x           = x,
  y           = y,
  node_id     = node_id,
  tree        = rj_tree,
  max_iter    = 5, alpha = 0.05,
  nsim        = 999, seed = 2016,
  max_pop_pct = 0.50, n_cores = 4L
)
print(seq_rj)


## ---- 5. Polygons + cluster membership ----
## RJ municipal boundaries ship with the package (IBGE Malhas Municipais),
## so no download is needed -- mirrors chicago_map / london_boroughs_map.
data(rj_map)
mun <- rj_map        # columns: ibge_code (6-digit), NAME, geometry

region_info <- unique(rj_mortality[, c("region_id", "ibge_code", "name")])

cr1    <- merge(get_cluster_regions(result_rj, n_clusters = 1, overlap = FALSE),
                region_info, by = "region_id")
cr2    <- merge(get_cluster_regions(result_rj, n_clusters = 2, overlap = TRUE),
                region_info, by = "region_id")
cr_seq <- merge(get_cluster_regions(seq_rj, overlap = TRUE),
                region_info, by = "region_id")


## ---- 6. Plot 1: Most likely cluster ----
mun1 <- merge(mun, cr1, by = "ibge_code", all.x = TRUE)
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
mun2 <- merge(mun, cr2, by = "ibge_code")
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


## ---- 8. Plot 3: Sequential scan clusters ----
n_iter <- seq_rj$n_iter
if (n_iter > 0) {
  # Cross-join the map polygons with the panels of the sequential scan,
  # so that every municipality appears in every panel (those that are
  # outside the dataset of 89 municipalities show up as NA-coloured).
  panels <- unique(cr_seq$panel)
  cr_seq_keys <- unique(cr_seq[, c("ibge_code", "cluster", "node_id",
                                    "llr", "pvalue", "panel")])
  mun_panels <- merge(
    do.call(rbind,
            lapply(panels, function(p) cbind(mun, panel = p))),
    cr_seq_keys,
    by = c("ibge_code", "panel"),
    all.x = TRUE
  )

  palette_seq <- c("#C44E52", "#4C72B0", "#55A868", "#8172B2", "#CCB974")[
    seq_len(n_iter)
  ]
  names(palette_seq) <- as.character(seq_len(n_iter))

  n_sig <- sum(seq_rj$clusters$significant, na.rm = TRUE)
  p_seq <- ggplot(mun_panels) +
    geom_sf(aes(fill = factor(cluster)), color = "gray40",
            linewidth = 0.12, alpha = 0.75) +
    scale_fill_manual(values = palette_seq, na.value = "gray95",
                      name = "Iteration", na.translate = FALSE) +
    facet_wrap(~ panel, nrow = 1) +
    labs(title    = paste0("Sequential Scan (Zhang et al. 2010): ",
                            n_iter, " iterations, ",
                            n_sig, " significant"),
         subtitle = "Each panel = one iteration (regions removed before next).") +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(color = "gray40"),
          strip.text    = element_text(face = "bold", size = 9),
          axis.title    = element_blank(),
          axis.text     = element_text(color = "gray50", size = 7),
          legend.position = "none")

  ggsave("rj_clusters_sequential.png", p_seq,
         width = max(8, 4 * n_iter), height = 6, dpi = 300)
  cat("Saved: rj_clusters_sequential.png (", n_iter, " iterations)\n")
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
cat("Sequential scan (Zhang et al. 2010): ", n_iter, " iterations, ",
    sum(seq_rj$clusters$significant, na.rm = TRUE),
    " significant\n", sep = "")
cat("\nDone!\n")
