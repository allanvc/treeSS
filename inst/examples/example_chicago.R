##=============================================================================
## Example 2: Crime in Chicago, USA (2023)
##
## Two analyses:
##   1. PAPER-FAITHFUL multi-cluster procedure (Cancado et al. 2025, Sec. 5.1.1):
##      single scan + filter_clusters / get_cluster_regions(n_clusters = N).
##   2. (Optional) iterative_scan with Holm-Bonferroni - NOT in paper.
##
## Population denominator: ACS 2020 5-year residential population per
## community area (replaces the bundled compositional 'population' column).
##=============================================================================

options(bitmapType = "cairo")

library(treeSS)
library(ggplot2)
library(sf)


## ---- 1. Load data ----
data(chicago_crimes)
data(chicago_tree)
data(chicago_map)

cat("=== Chicago Crime 2023 ===\n")
cat("Rows (long):           ", nrow(chicago_crimes), "\n")
cat("Unique community areas:", length(unique(chicago_crimes$region_id)), "\n")
cat("Tree nodes:            ", nrow(chicago_tree), "\n")
cat("Total incidents:       ", sum(chicago_crimes$cases), "\n")

crime_desc <- function(node_id) {
  parts <- strsplit(node_id, " \\| ")[[1]]
  if (length(parts) == 3) return(paste0(parts[2], " (", parts[3], ")"))
  if (length(parts) == 2) return(parts[2])
  node_id
}


## ---- 2. Choose denominator: residents (not total incidents) ----
##
## chicago_crimes ships with two denominator columns:
##   * population:      total incidents per area (compositional - useful
##                       for "which crime types over-occur where", but
##                       not a population at risk)
##   * pop_residential: ACS 2020 5-year residential population. This is
##                       the appropriate denominator for an incidence-rate
##                       analysis (incidents per resident), analogous to
##                       deaths/live_births in the paper's RJ application.
##
## We use pop_residential below.
cat("Total Chicago residents:",
    sum(unique(chicago_crimes[, c("area_number", "pop_residential")])$pop_residential),
    "\n")


## ---- 3. Run the tree-spatial scan ----
cat("\nRunning tree-spatial scan (nsim=999, n_cores=4)...\n")
system.time({
  result_chi <- treespatial_scan(
    cases       = chicago_crimes$cases,
    population  = chicago_crimes$pop_residential,
    region_id   = chicago_crimes$region_id,
    x           = chicago_crimes$x,
    y           = chicago_crimes$y,
    node_id     = chicago_crimes$node_id,
    tree        = chicago_tree,
    max_pop_pct = 0.25,
    nsim        = 999, seed = 2023,
    n_cores     = 4L
  )
})
print(result_chi)


## ---- 4. Paper-faithful: distinct top clusters via filter_clusters ----
cat("\n=== Distinct top clusters (paper Sec. 5.1.1) ===\n")
fc <- filter_clusters(result_chi)
print(head(fc[, c("node_id", "n_regions", "cases", "expected", "llr", "pvalue")],
            5))


## ---- 5. (Optional) Iterative scan with Holm-Bonferroni ----
cat("\n=== Iterative scan (extension, with Holm-Bonferroni) ===\n")
iter_chi <- iterative_scan(
  cases       = chicago_crimes$cases,
  population  = chicago_crimes$pop_residential,
  region_id   = chicago_crimes$region_id,
  x           = chicago_crimes$x,
  y           = chicago_crimes$y,
  node_id     = chicago_crimes$node_id,
  tree        = chicago_tree,
  max_iter    = 5, alpha = 0.05,
  nsim        = 999, seed = 2023,
  max_pop_pct = 0.25, n_cores = 4L
)
print(iter_chi)


## ---- 6. Cluster membership + descriptive lookup ----
region_info <- unique(chicago_crimes[, c("region_id", "area_number", "name")])

cr1   <- merge(get_cluster_regions(result_chi, n_clusters = 1, overlap = FALSE),
               region_info, by = "region_id")
cr2   <- merge(get_cluster_regions(result_chi, n_clusters = 2, overlap = TRUE),
               region_info, by = "region_id")
cr_it <- merge(get_cluster_regions(iter_chi, overlap = TRUE),
               region_info, by = "region_id")


## ---- 7. Plot 1: Most likely cluster ----
chi1 <- merge(chicago_map, cr1, by.x = "AREA_NUM", by.y = "area_number",
              all.x = TRUE)
mlc <- result_chi$most_likely_cluster

p1 <- ggplot(chi1) +
  geom_sf(aes(fill = factor(cluster)), color = "gray40",
          linewidth = 0.15, alpha = 0.75) +
  scale_fill_manual(
    values   = c("1" = "#C44E52"),
    na.value = "gray95",
    labels   = c("1" = paste0(crime_desc(mlc$node_id), " (",
                              length(mlc$region_ids), " areas)")),
    name     = "Cluster"
  ) +
  labs(
    title    = "Tree-Spatial Scan: most likely cluster",
    subtitle = paste0(crime_desc(mlc$node_id),
                      " | LR=", round(mlc$llr, 1),
                      ", p=", format.pval(result_chi$pvalue, digits = 3),
                      " | ", mlc$cases, " incidents (expected ",
                      round(mlc$expected, 1), ")")
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray40", size = 10),
        axis.title    = element_blank(),
        axis.text     = element_text(color = "gray50", size = 8))

ggsave("chicago_cluster_mlc.png", p1, width = 9, height = 9, dpi = 300)
cat("Saved: chicago_cluster_mlc.png\n")


## ---- 8. Plot 2: Top-2 distinct clusters (paper Sec. 5.1.1) ----
chi2 <- merge(chicago_map, cr2, by.x = "AREA_NUM", by.y = "area_number")
palette2 <- c("1" = "#C44E52", "2" = "#4C72B0")

p2 <- ggplot(chi2) +
  geom_sf(aes(fill = factor(cluster)), color = "gray40",
          linewidth = 0.12, alpha = 0.75) +
  scale_fill_manual(values = palette2, na.value = "gray95",
                    name = "Cluster", na.translate = FALSE) +
  facet_wrap(~ panel, nrow = 1) +
  labs(title    = "Tree-Spatial Scan: Top 2 Distinct Clusters",
       subtitle = "Cancado et al. (2025), Sec. 5.1.1 - denominator: residents") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray40"),
        strip.text    = element_text(face = "bold", size = 9),
        axis.title    = element_blank(),
        axis.text     = element_text(color = "gray50", size = 7),
        legend.position = "none")

ggsave("chicago_clusters_top2.png", p2, width = 14, height = 6, dpi = 300)
cat("Saved: chicago_clusters_top2.png\n")


## ---- 9. Plot 3 (optional): Iterative scan clusters ----
n_iter <- iter_chi$n_iter
if (n_iter > 0) {
  chi_it <- merge(chicago_map, cr_it,
                  by.x = "AREA_NUM", by.y = "area_number", all.x = TRUE)
  palette_it <- c("#C44E52", "#4C72B0", "#55A868", "#8172B2", "#CCB974")[
    seq_len(n_iter)
  ]
  names(palette_it) <- as.character(seq_len(n_iter))

  n_sig <- sum(iter_chi$clusters$significant, na.rm = TRUE)
  p_it <- ggplot(chi_it) +
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

  ggsave("chicago_clusters_iterative.png", p_it,
         width = max(8, 4 * n_iter), height = 6, dpi = 300)
  cat("Saved: chicago_clusters_iterative.png (", n_iter, " iterations)\n")
}


## ---- 10. Summary ----
cat("\n=== SUMMARY ===\n")
cat("Most likely cluster:\n")
cat("  Node:", mlc$node_id, "\n")
cat("  Description:", crime_desc(mlc$node_id), "\n")
cat("  Areas:", length(mlc$region_ids), "\n")
cat("  Incidents:", mlc$cases, "(expected", round(mlc$expected, 1), ")\n")
cat("  RR:", round(mlc$rr, 2),
    "  LR:", round(mlc$llr, 2),
    "  p-value:", result_chi$pvalue, "\n")
cat("\nDistinct clusters via filter_clusters:", nrow(fc), "\n")
cat("Iterative scan: ", n_iter, " iterations, ",
    sum(iter_chi$clusters$significant, na.rm = TRUE),
    " significant after Holm-Bonferroni\n", sep = "")
cat("\nDone!\n")
