##=============================================================================
## Example 4: General Mortality in Florida, USA (2016)
##
## Two analyses:
##   1. PAPER-FAITHFUL multi-cluster procedure (Cancado et al. 2025, Sec. 5.1.1):
##      single scan + filter_clusters / get_cluster_regions(n_clusters = N).
##   2. (Optional) iterative_scan with Holm-Bonferroni - NOT in paper.
##
## Pedagogical workflow:
##   - Load raw long-format data (fl_deaths)
##   - Build ICD-10 tree directly from the codes that actually appear
##   - Download county polygons + centroids from tigris
##   - Build parallel vectors (cases, population, region_id, x, y, node_id)
##   - Run scan + iterative scan with Holm-Bonferroni
##
## Data: CDC WONDER Compressed Mortality File 1999-2016
##=============================================================================

options(bitmapType = "cairo")

library(treeSS)
library(ggplot2)
library(sf)


## ---- 1. Load raw data ----
data(fl_deaths)

cat("=== Florida Mortality 2016 ===\n")
cat("Rows (raw):  ", nrow(fl_deaths), "\n")
cat("Counties:    ", length(unique(fl_deaths$county_fips)), "\n")
cat("ICD-10 codes:", length(unique(fl_deaths$icd10_code)), "\n")
cat("Total deaths:", sum(fl_deaths$deaths), "\n")

icd_lookup <- unique(fl_deaths[, c("icd10_code", "icd10_desc")])
icd_desc <- function(code) {
  d <- icd_lookup$icd10_desc[icd_lookup$icd10_code == code]
  if (length(d) == 0) return(code)
  d[1]
}


## ---- 2. Build ICD-10 tree directly from fl_deaths ----
##
## fl_deaths uses two kinds of codes:
##   * 3-char (e.g. "C56", "R54") - leaf at the chapter level
##   * 5-char with a dot (e.g. "A04.7") - leaf under a 3-char internal parent

cat("\nBuilding ICD-10 tree from fl_deaths codes...\n")
icd_codes <- sort(unique(fl_deaths$icd10_code))

parsed <- data.frame(
  code = icd_codes,
  three = ifelse(grepl("\\.", icd_codes),
                  sub("\\..*", "", icd_codes),
                  icd_codes),
  is_three_leaf = !grepl("\\.", icd_codes) & nchar(icd_codes) == 3,
  stringsAsFactors = FALSE
)

icd_chapter_of <- function(three) {
  L <- substr(three, 1, 1)
  num <- suppressWarnings(as.integer(substr(three, 2, 3)))
  switch(
    L,
    A = "I_AB", B = "I_AB",
    C = "II_C_D",
    D = if (!is.na(num) && num <= 48) "II_C_D" else "III_D",
    E = "IV_E",
    F = "V_F",
    G = "VI_G",
    H = if (!is.na(num) && num <= 59) "VII_H" else "VIII_H",
    I = "IX_I",
    J = "X_J",
    K = "XI_K",
    L = "XII_L",
    M = "XIII_M",
    N = "XIV_N",
    O = "XV_O",
    P = "XVI_P",
    Q = "XVII_Q",
    R = "XVIII_R",
    S = "XIX_S_T", T = "XIX_S_T",
    V = "XX_V_Y", W = "XX_V_Y", X = "XX_V_Y", Y = "XX_V_Y",
    "Other"
  )
}
parsed$chapter <- sapply(parsed$three, icd_chapter_of)

fl_tree <- unique(rbind(
  data.frame(node_id = "Root", parent_id = NA_character_,
             stringsAsFactors = FALSE),
  data.frame(node_id = unique(parsed$chapter), parent_id = "Root",
             stringsAsFactors = FALSE),
  unique(data.frame(node_id = parsed$three, parent_id = parsed$chapter,
                    stringsAsFactors = FALSE)),
  data.frame(
    node_id   = parsed$code[!parsed$is_three_leaf],
    parent_id = parsed$three[!parsed$is_three_leaf],
    stringsAsFactors = FALSE
  )
))

parents <- unique(fl_tree$parent_id[!is.na(fl_tree$parent_id)])
leaves  <- setdiff(fl_tree$node_id, parents)
stopifnot(length(setdiff(icd_codes, leaves)) == 0)
cat("  Total nodes:", nrow(fl_tree),
    "| Chapters:", length(unique(parsed$chapter)),
    "| Leaves:", length(leaves), "\n")


## ---- 3. Download county polygons + centroids ----
if (!requireNamespace("tigris", quietly = TRUE)) {
  stop("Install 'tigris': install.packages('tigris')")
}
library(tigris)
cat("\nDownloading FL county boundaries...\n")
fl_map <- counties(state = "FL", cb = TRUE, year = 2016)
fl_map <- st_transform(fl_map, 4326)
centroids <- st_centroid(st_geometry(fl_map))
coords <- st_coordinates(centroids)
fl_map$x <- coords[, "X"]
fl_map$y <- coords[, "Y"]


## ---- 4. Build parallel vectors ----
fl_pop <- aggregate(population ~ county_fips, data = fl_deaths, FUN = max)
xy <- data.frame(county_fips = fl_map$GEOID, x = fl_map$x, y = fl_map$y,
                  stringsAsFactors = FALSE)
dat <- merge(fl_deaths[, c("county_fips", "icd10_code", "deaths")],
             fl_pop, by = "county_fips")
dat <- merge(dat, xy, by = "county_fips")
dat$region_id <- as.integer(factor(dat$county_fips,
                                    levels = sort(unique(dat$county_fips))))
dat <- dat[dat$deaths > 0, ]

region_info <- unique(dat[, c("region_id", "county_fips")])
cat("Long rows:", nrow(dat), "\n")


## ---- 5. Run the tree-spatial scan ----
cat("\nRunning tree-spatial scan (nsim=999, n_cores=4, max_pop_pct=0.05)...\n")
system.time({
  result_fl <- treespatial_scan(
    cases       = dat$deaths,
    population  = dat$population,
    region_id   = dat$region_id,
    x           = dat$x,
    y           = dat$y,
    node_id     = dat$icd10_code,
    tree        = fl_tree,
    max_pop_pct = 0.05,
    nsim        = 999, seed = 2016,
    n_cores     = 4L
  )
})
print(result_fl)


## ---- 6. Paper-faithful: distinct top clusters via filter_clusters ----
cat("\n=== Distinct top clusters (paper Sec. 5.1.1) ===\n")
fc <- filter_clusters(result_fl)
print(head(fc[, c("node_id", "n_regions", "cases", "expected", "llr", "pvalue")],
            5))


## ---- 7. (Optional) Iterative scan with Holm-Bonferroni ----
cat("\n=== Iterative scan (extension, with Holm-Bonferroni) ===\n")
iter_fl <- iterative_scan(
  cases       = dat$deaths,
  population  = dat$population,
  region_id   = dat$region_id,
  x           = dat$x,
  y           = dat$y,
  node_id     = dat$icd10_code,
  tree        = fl_tree,
  max_iter    = 5, alpha = 0.05,
  nsim        = 999, seed = 2016,
  max_pop_pct = 0.05, n_cores = 4L
)
print(iter_fl)


## ---- 8. Cluster membership ----
cr1   <- merge(get_cluster_regions(result_fl, n_clusters = 1, overlap = FALSE),
               region_info, by = "region_id")
cr2   <- merge(get_cluster_regions(result_fl, n_clusters = 2, overlap = TRUE),
               region_info, by = "region_id")
cr_it <- merge(get_cluster_regions(iter_fl, overlap = TRUE),
               region_info, by = "region_id")


## ---- 9. Plot 1: Most likely cluster ----
fl1 <- merge(fl_map, cr1, by.x = "GEOID", by.y = "county_fips", all.x = TRUE)
mlc <- result_fl$most_likely_cluster

p1 <- ggplot(fl1) +
  geom_sf(aes(fill = factor(cluster)), color = "gray40",
          linewidth = 0.15, alpha = 0.75) +
  scale_fill_manual(
    values   = c("1" = "#C44E52"),
    na.value = "gray95",
    labels   = c("1" = paste0(mlc$node_id, " (",
                              length(mlc$region_ids), " counties)")),
    name     = "Cluster"
  ) +
  labs(
    title    = paste0("Tree-Spatial Scan: most likely cluster (", mlc$node_id, ")"),
    subtitle = paste0(icd_desc(mlc$node_id),
                      " | LR=", round(mlc$llr, 1),
                      ", p=", format.pval(result_fl$pvalue, digits = 3),
                      " | ", mlc$cases, " deaths (expected ",
                      round(mlc$expected, 1), ")")
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray40", size = 10),
        axis.title    = element_blank(),
        axis.text     = element_text(color = "gray50", size = 8))

ggsave("fl_cluster_mlc.png", p1, width = 9, height = 8, dpi = 300)
cat("Saved: fl_cluster_mlc.png\n")


## ---- 10. Plot 2: Top-2 distinct clusters (paper Sec. 5.1.1) ----
fl2 <- merge(fl_map, cr2, by.x = "GEOID", by.y = "county_fips")
palette2 <- c("1" = "#C44E52", "2" = "#4C72B0")

p2 <- ggplot(fl2) +
  geom_sf(aes(fill = factor(cluster)), color = "gray40",
          linewidth = 0.12, alpha = 0.75) +
  scale_fill_manual(values = palette2, na.value = "gray95",
                    name = "Cluster", na.translate = FALSE) +
  facet_wrap(~ panel, nrow = 1) +
  labs(title    = "Tree-Spatial Scan: Top 2 Distinct Clusters",
       subtitle = "Cancado et al. (2025), Sec. 5.1.1 - Florida mortality, 2016") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray40"),
        strip.text    = element_text(face = "bold", size = 10),
        axis.title    = element_blank(),
        axis.text     = element_text(color = "gray50", size = 7),
        legend.position = "none")

ggsave("fl_clusters_top2.png", p2, width = 14, height = 6, dpi = 300)
cat("Saved: fl_clusters_top2.png\n")


## ---- 11. Plot 3 (optional): Iterative scan clusters ----
n_iter <- iter_fl$n_iter
if (n_iter > 0) {
  fl_it <- merge(fl_map, cr_it, by.x = "GEOID", by.y = "county_fips", all.x = TRUE)
  palette_it <- c("#C44E52", "#4C72B0", "#55A868", "#8172B2", "#CCB974")[
    seq_len(n_iter)
  ]
  names(palette_it) <- as.character(seq_len(n_iter))

  n_sig <- sum(iter_fl$clusters$significant, na.rm = TRUE)
  p_it <- ggplot(fl_it) +
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

  ggsave("fl_clusters_iterative.png", p_it,
         width = max(8, 4 * n_iter), height = 6, dpi = 300)
  cat("Saved: fl_clusters_iterative.png (", n_iter, " iterations)\n")
}


## ---- 12. Summary ----
cat("\n=== SUMMARY ===\n")
cat("Most likely cluster:\n")
cat("  Node:", mlc$node_id, "(", icd_desc(mlc$node_id), ")\n")
cat("  Counties:", length(mlc$region_ids), "\n")
cat("  Deaths:", mlc$cases, "(expected", round(mlc$expected, 1), ")\n")
cat("  RR:", round(mlc$rr, 2),
    "  LR:", round(mlc$llr, 2),
    "  p-value:", result_fl$pvalue, "\n")
cat("\nDistinct clusters via filter_clusters:", nrow(fc), "\n")
cat("Iterative scan: ", n_iter, " iterations, ",
    sum(iter_fl$clusters$significant, na.rm = TRUE),
    " significant after Holm-Bonferroni\n", sep = "")
cat("\nDone!\n")
