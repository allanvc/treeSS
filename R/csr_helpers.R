# CSR format conversion helpers for the C++ backend

#' @keywords internal
.zones_to_csr <- function(zones) {
  # Convert zones list to CSR (compressed sparse row) format for C++
  # Returns: zone_region_idx (0-based), zone_ptr, zone_pop
  n_zones <- length(zones)

  lens <- vapply(zones, function(z) length(z$region_idx), integer(1))
  zone_ptr <- c(0L, cumsum(lens))

  if (n_zones > 0L && sum(lens) > 0L) {
    all_idx <- unlist(lapply(zones, `[[`, "region_idx"),
                      use.names = FALSE) - 1L  # 0-based for C++
  } else {
    all_idx <- integer(0)
  }

  zone_pop <- vapply(zones, function(z) as.numeric(z$population), numeric(1))

  list(
    zone_region_idx = as.integer(all_idx),
    zone_ptr = as.integer(zone_ptr),
    zone_pop = as.numeric(zone_pop)
  )
}

#' @keywords internal
.tree_to_csr_children <- function(tree) {
  # Convert tree children to CSR format for C++
  # Returns: children_idx (0-based), children_ptr
  all_nodes <- tree$node_id
  n_nodes <- length(all_nodes)

  # Row position of each node's parent within all_nodes (NA for the root and
  # for any parent_id not present in node_id, matching the previous
  # behaviour of dropping such rows).
  parent_pos <- match(tree$parent_id, all_nodes)
  child_rows <- which(!is.na(parent_pos))

  if (length(child_rows) > 0L) {
    # Group children by parent position; within a parent, keep increasing
    # row order (same order the previous per-node which() scan produced).
    ord <- order(parent_pos[child_rows], child_rows)
    children_sorted <- child_rows[ord]
    counts <- tabulate(parent_pos[child_rows], nbins = n_nodes)
  } else {
    children_sorted <- integer(0)
    counts <- integer(n_nodes)
  }

  children_ptr <- c(0L, cumsum(counts))

  list(
    children_idx = as.integer(children_sorted - 1L),  # 0-based for C++
    children_ptr = as.integer(children_ptr)
  )
}
