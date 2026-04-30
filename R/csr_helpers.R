# CSR format conversion helpers for the C++ backend

#' @keywords internal
.zones_to_csr <- function(zones) {
  # Convert zones list to CSR (compressed sparse row) format for C++
  # Returns: zone_region_idx (0-based), zone_ptr, zone_pop
  n_zones <- length(zones)
  all_idx <- integer(0)
  zone_ptr <- integer(n_zones + 1)
  zone_pop <- numeric(n_zones)

  ptr <- 0L
  for (zi in seq_len(n_zones)) {
    zone_ptr[zi] <- ptr
    idx <- zones[[zi]]$region_idx - 1L  # 0-based for C++
    all_idx <- c(all_idx, idx)
    ptr <- ptr + length(idx)
    zone_pop[zi] <- zones[[zi]]$population
  }
  zone_ptr[n_zones + 1] <- ptr

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
  all_idx <- integer(0)
  children_ptr <- integer(n_nodes + 1)

  ptr <- 0L
  for (i in seq_len(n_nodes)) {
    children_ptr[i] <- ptr
    ch <- which(!is.na(tree$parent_id) & tree$parent_id == all_nodes[i])
    if (length(ch) > 0) {
      all_idx <- c(all_idx, ch - 1L)  # 0-based for C++
      ptr <- ptr + length(ch)
    }
  }
  children_ptr[n_nodes + 1] <- ptr

  list(
    children_idx = as.integer(all_idx),
    children_ptr = as.integer(children_ptr)
  )
}
