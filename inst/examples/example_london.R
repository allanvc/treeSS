##=============================================================================
## Example 3: Road Traffic Collisions in London, UK (2022)
##
## Two analyses:
##   1. PAPER-FAITHFUL multi-cluster procedure (Cancado et al. 2025, Sec. 5.1.1):
##      single scan + filter_clusters / get_cluster_regions(n_clusters = N).
##   2. (Optional) iterative_scan with Holm-Bonferroni - NOT in paper.
##
## Maps: leaflet (interactive, HTML)
## Polygons: london_boroughs_map (included in the package)
##=============================================================================

library(treeSS)
library(leaflet)
library(sf)
library(htmlwidgets)
library(htmltools)


## ---- 1. Load data ----
data(london_collisions)
data(london_tree)
data(london_boroughs_map)

cat("=== London Road Collisions 2022 ===\n")
cat("Rows (long):     ", nrow(london_collisions), "\n")
cat("Unique boroughs: ", length(unique(london_collisions$region_id)), "\n")
cat("Tree nodes:      ", nrow(london_tree), "\n")
cat("Total collisions:", sum(london_collisions$cases), "\n")


## ---- 2. Run the tree-spatial scan ----
cat("\nRunning tree-spatial scan (nsim=999, n_cores=4)...\n")
system.time({
  result_ldn <- treespatial_scan(
    cases       = london_collisions$cases,
    population  = london_collisions$population,
    region_id   = london_collisions$region_id,
    x           = london_collisions$x,
    y           = london_collisions$y,
    node_id     = london_collisions$node_id,
    tree        = london_tree,
    max_pop_pct = 0.25,
    nsim        = 999, seed = 42,
    n_cores     = 4L
  )
})
print(result_ldn)


## ---- 3. Paper-faithful: distinct top clusters via filter_clusters ----
cat("\n=== Distinct top clusters (paper Sec. 5.1.1) ===\n")
fc <- filter_clusters(result_ldn)
print(head(fc[, c("node_id", "n_regions", "cases", "expected", "llr", "pvalue")],
            5))


## ---- 4. (Optional) Iterative scan with Holm-Bonferroni ----
cat("\n=== Iterative scan (extension, with Holm-Bonferroni) ===\n")
iter_ldn <- iterative_scan(
  cases       = london_collisions$cases,
  population  = london_collisions$population,
  region_id   = london_collisions$region_id,
  x           = london_collisions$x,
  y           = london_collisions$y,
  node_id     = london_collisions$node_id,
  tree        = london_tree,
  max_iter    = 5, alpha = 0.05,
  nsim        = 999, seed = 42,
  max_pop_pct = 0.25, n_cores = 4L
)
print(iter_ldn)


## ---- 5. Cluster membership + descriptive lookup ----
region_info <- unique(london_collisions[, c("region_id", "borough")])

cr1   <- merge(get_cluster_regions(result_ldn, n_clusters = 1, overlap = FALSE),
               region_info, by = "region_id")
cr2   <- merge(get_cluster_regions(result_ldn, n_clusters = 2, overlap = TRUE),
               region_info, by = "region_id")
cr_it <- merge(get_cluster_regions(iter_ldn, overlap = TRUE),
               region_info, by = "region_id")

JOIN_X <- "NAME"     # column in london_boroughs_map
JOIN_Y <- "borough"  # column in cluster results


## ---- 6. Map 1: Most likely cluster (leaflet) ----
ldn1 <- merge(london_boroughs_map, cr1,
              by.x = JOIN_X, by.y = JOIN_Y, all.x = TRUE)
mlc <- result_ldn$most_likely_cluster

ldn1$label <- ifelse(
  is.na(ldn1$cluster),
  paste0("<b>", ldn1$NAME, "</b><br>Not in cluster"),
  paste0("<b>", ldn1$NAME, "</b><br>",
         "<span style='color:#C44E52'><b>", ldn1$node_id, "</b></span><br>",
         "LR: ", round(ldn1$llr, 1),
         " | p = ", round(ldn1$pvalue, 3))
)

map1 <- leaflet(ldn1) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(
    fillColor   = ~ifelse(is.na(cluster), "#E8E8E8", "#C44E52"),
    fillOpacity = ~ifelse(is.na(cluster), 0.3, 0.65),
    color       = "gray40", weight = 1.5,
    label       = ~lapply(label, HTML)
  ) %>%
  setView(lng = -0.12, lat = 51.5, zoom = 10) %>%
  addLegend(
    position = "bottomright",
    colors   = c("#C44E52", "#E8E8E8"),
    labels   = c(paste0(mlc$node_id, " (",
                        length(mlc$region_ids), " boroughs)"),
                 "Not in cluster"),
    title    = paste0("LR=", round(mlc$llr, 1),
                      ", p=", format.pval(result_ldn$pvalue, digits = 3)),
    opacity  = 0.8
  )

saveWidget(map1, "london_cluster_mlc.html", selfcontained = TRUE)
cat("Saved: london_cluster_mlc.html\n")


## ---- 7. Map 2: Top-2 distinct clusters (paper Sec. 5.1.1) ----
# A borough may appear in multiple cluster panels in overlap mode; for the
# leaflet map we keep the lowest-ranked cluster (cluster 1 wins over 2).
cr2_first <- cr2[order(cr2$cluster), ]
cr2_first <- cr2_first[!duplicated(cr2_first$borough), ]

ldn2 <- merge(london_boroughs_map, cr2_first,
              by.x = JOIN_X, by.y = JOIN_Y, all.x = TRUE)

palette2 <- c("1" = "#C44E52", "2" = "#4C72B0")
top2 <- filter_clusters(result_ldn)
top2 <- head(top2, 2)

ldn2$fillc <- ifelse(is.na(ldn2$cluster), "#E8E8E8",
                     palette2[as.character(ldn2$cluster)])
ldn2$label <- ifelse(
  is.na(ldn2$cluster),
  paste0("<b>", ldn2$NAME, "</b><br>Not in any top-2 cluster"),
  paste0("<b>", ldn2$NAME, "</b><br>",
         "Cluster ", ldn2$cluster, ": ",
         "<span style='color:", palette2[as.character(ldn2$cluster)],
         "'><b>", ldn2$node_id, "</b></span><br>",
         "LR: ", round(ldn2$llr, 1),
         " | p = ", round(ldn2$pvalue, 3))
)

map2 <- leaflet(ldn2) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(
    fillColor   = ~fillc,
    fillOpacity = ~ifelse(is.na(cluster), 0.3, 0.65),
    color       = "gray40", weight = 1.5,
    label       = ~lapply(label, HTML)
  ) %>%
  setView(lng = -0.12, lat = 51.5, zoom = 10) %>%
  addLegend(
    position = "bottomright",
    colors   = c(palette2, "#E8E8E8"),
    labels   = c(paste0("Cluster ", seq_along(top2$node_id), ": ",
                        top2$node_id),
                 "Other"),
    title    = "Top 2 distinct (paper Sec. 5.1.1)",
    opacity  = 0.8
  )

saveWidget(map2, "london_clusters_top2.html", selfcontained = TRUE)
cat("Saved: london_clusters_top2.html\n")


## ---- 8. Map 3 (optional): Iterative scan clusters ----
n_iter <- iter_ldn$n_iter
if (n_iter > 0) {
  cr_it_first <- cr_it[order(cr_it$cluster), ]
  cr_it_first <- cr_it_first[!duplicated(cr_it_first$borough), ]

  ldn_it <- merge(london_boroughs_map, cr_it_first,
                  by.x = JOIN_X, by.y = JOIN_Y, all.x = TRUE)

  palette_it <- c("#C44E52", "#4C72B0", "#55A868", "#8172B2", "#CCB974")[
    seq_len(n_iter)
  ]
  names(palette_it) <- as.character(seq_len(n_iter))

  n_sig <- sum(iter_ldn$clusters$significant, na.rm = TRUE)
  ldn_it$fillc <- ifelse(is.na(ldn_it$cluster), "#E8E8E8",
                          palette_it[as.character(ldn_it$cluster)])
  ldn_it$label <- ifelse(
    is.na(ldn_it$cluster),
    paste0("<b>", ldn_it$NAME, "</b><br>Not detected"),
    paste0("<b>", ldn_it$NAME, "</b><br>",
           "Iteration ", ldn_it$cluster, ": ",
           "<span style='color:", palette_it[as.character(ldn_it$cluster)],
           "'><b>", ldn_it$node_id, "</b></span><br>",
           "LR: ", round(ldn_it$llr, 1),
           " | adj p = ", round(ldn_it$pvalue_adjusted, 3))
  )

  map_it <- leaflet(ldn_it) %>%
    addProviderTiles(providers$CartoDB.Positron) %>%
    addPolygons(
      fillColor   = ~fillc,
      fillOpacity = ~ifelse(is.na(cluster), 0.3, 0.65),
      color       = "gray40", weight = 1.5,
      label       = ~lapply(label, HTML)
    ) %>%
    setView(lng = -0.12, lat = 51.5, zoom = 10) %>%
    addLegend(
      position = "bottomright",
      colors   = c(palette_it, "#E8E8E8"),
      labels   = c(paste0("Iteration ", seq_len(n_iter)), "Other"),
      title    = paste0("Iterative scan (extension)<br>",
                        n_sig, " of ", n_iter,
                        " significant (Holm-Bonferroni)"),
      opacity  = 0.8
    )

  saveWidget(map_it, "london_clusters_iterative.html", selfcontained = TRUE)
  cat("Saved: london_clusters_iterative.html\n")
}


## ---- 9. Summary ----
cat("\n=== SUMMARY ===\n")
cat("Most likely cluster:\n")
cat("  Node:", mlc$node_id, "\n")
cat("  Boroughs:", length(mlc$region_ids), "\n")
cat("  Collisions:", mlc$cases, "(expected", round(mlc$expected, 1), ")\n")
cat("  RR:", round(mlc$rr, 2),
    "  LR:", round(mlc$llr, 2),
    "  p-value:", result_ldn$pvalue, "\n")
cat("\nDistinct clusters via filter_clusters:", nrow(fc), "\n")
cat("Iterative scan: ", n_iter, " iterations, ",
    sum(iter_ldn$clusters$significant, na.rm = TRUE),
    " significant after Holm-Bonferroni\n", sep = "")
cat("\nDone!\n")
