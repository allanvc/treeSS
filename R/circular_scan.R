#' Kulldorff's Circular Spatial Scan Statistic
#'
#' Performs Kulldorff's circular spatial scan statistic for detecting
#' spatial clusters. Inputs are passed as parallel vectors with one entry
#' per region (cases must already be aggregated to the region level).
#'
#' @param cases Numeric vector of length \eqn{n} (one entry per region):
#'   total cases in each region.
#' @param population Numeric vector of length \eqn{n}: population (or
#'   denominator) of each region.
#' @param region_id Vector of region identifiers, length \eqn{n}.
#' @param x,y Numeric vectors of region centroid coordinates, length \eqn{n}.
#' @param max_pop_pct Numeric. Default \code{0.5}.
#' @param nsim Integer. Number of MC simulations. Default \code{999}.
#' @param alpha Numeric. Significance level. Default \code{0.05}.
#' @param n_secondary Integer. Default \code{1000}.
#' @param model Character. \code{"poisson"} or \code{"binomial"}.
#' @param seed Integer or \code{NULL}.
#' @param n_cores Integer. OpenMP threads.
#'
#' @return An object of class \code{"circular_scan"}.
#'
#' @references
#' Kulldorff, M. (1997). A spatial scan statistic. \emph{Communications in
#' Statistics - Theory and Methods}, 26(6), 1481-1496.
#'
#' @seealso \code{\link{filter_clusters}}, \code{\link{tree_scan}},
#'   \code{\link{treespatial_scan}}, \code{\link{get_cluster_regions}},
#'   \code{\link{iterative_scan}}
#'
#' @export
#' @examples
#' set.seed(42)
#' n <- 20
#' cases <- rpois(n, lambda = 10)
#' cases[1:5] <- rpois(5, lambda = 30)
#'
#' result <- circular_scan(
#'   cases       = cases,
#'   population  = rep(1000, n),
#'   region_id   = 1:n,
#'   x           = runif(n, 0, 10),
#'   y           = runif(n, 0, 10),
#'   nsim        = 99
#' )
#' print(result)
circular_scan <- function(cases, population, region_id, x, y,
                          max_pop_pct = 0.5, nsim = 999L,
                          alpha = 0.05, n_secondary = 1000L,
                          model = c("poisson", "binomial"),
                          seed = NULL, n_cores = 1L) {

  model <- match.arg(model)
  model_int <- if (model == "binomial") 1L else 0L

  regions <- .build_regions_circular(cases, population, region_id, x, y)

  .validate_regions(regions)

  if (!is.null(seed)) set.seed(seed)

  n <- nrow(regions)
  N <- sum(regions$population)
  C <- sum(cases)
  max_pop <- max_pop_pct * N

  zones <- build_zones(regions, max_pop = max_pop)

  if (length(zones) == 0) {
    stop("No valid candidate zones found.", call. = FALSE)
  }

  # C++ backend for Monte Carlo
  prob <- regions$population / N
  zone_csr <- .zones_to_csr(zones)

  cpp_result <- mc_spatial_cpp(
    as.numeric(cases),
    as.numeric(N),
    as.numeric(C),
    zone_csr$zone_region_idx,
    zone_csr$zone_ptr,
    zone_csr$zone_pop,
    as.numeric(prob),
    as.integer(nsim),
    model_int,
    as.integer(max(1L, n_cores))
  )

  obs_llr <- cpp_result$obs_max_llr
  best_z <- cpp_result$obs_best_z  # 1-based
  sim_llr <- cpp_result$sim_llr

  pvalue <- (sum(sim_llr >= obs_llr) + 1) / (nsim + 1)

  best_idx <- zones[[best_z]]$region_idx
  best_pop <- zones[[best_z]]$population
  best_cases <- sum(cases[best_idx])
  expected <- C * best_pop / N

  # --- Compute LLR for all zones (for secondary clusters) ---
  n_zones <- length(zones)
  zone_cases <- numeric(n_zones)
  zone_pop_vec <- numeric(n_zones)
  zone_centers <- integer(n_zones)
  zone_sizes <- integer(n_zones)

  for (z in seq_len(n_zones)) {
    idx <- zones[[z]]$region_idx
    zone_cases[z] <- sum(cases[idx])
    zone_pop_vec[z] <- zones[[z]]$population
    zone_centers[z] <- zones[[z]]$center
    zone_sizes[z] <- length(idx)
  }

  zone_expected <- C * zone_pop_vec / N
  zone_cases_bar <- C - zone_cases
  zone_expected_bar <- C - zone_expected
  zone_pop_bar <- N - zone_pop_vec

  if (model == "binomial") {
    rate_in  <- zone_cases / pmax(zone_pop_vec, 1e-300)
    rate_out <- zone_cases_bar / pmax(zone_pop_bar, 1e-300)
    mask <- rate_in > rate_out & zone_cases > 0
    p0 <- C / N

    t1  <- ifelse(mask & zone_cases > 0,
                  zone_cases * log(pmax(rate_in, 1e-300)), 0)
    t1b <- ifelse(mask & (zone_pop_vec - zone_cases) > 0,
                  (zone_pop_vec - zone_cases) * log(pmax(1 - rate_in, 1e-300)), 0)
    t2  <- ifelse(mask & zone_cases_bar > 0,
                  zone_cases_bar * log(pmax(rate_out, 1e-300)), 0)
    t2b <- ifelse(mask & (zone_pop_bar - zone_cases_bar) > 0,
                  (zone_pop_bar - zone_cases_bar) * log(pmax(1 - rate_out, 1e-300)), 0)
    t0  <- C * log(pmax(p0, 1e-300)) + (N - C) * log(pmax(1 - p0, 1e-300))
    zone_llr <- ifelse(mask, t1 + t1b + t2 + t2b - t0, 0)
    zone_llr[zone_llr < 0] <- 0
  } else {
    mask <- zone_cases > zone_expected & zone_cases > 0 & zone_expected > 0
    t1 <- ifelse(mask, zone_cases * log(pmax(zone_cases, 1e-300) /
                                         pmax(zone_expected, 1e-300)), 0)
    t2 <- ifelse(mask & zone_cases_bar > 0,
                 zone_cases_bar * log(pmax(zone_cases_bar, 1e-300) /
                                       pmax(zone_expected_bar, 1e-300)), 0)
    zone_llr <- ifelse(mask, t1 + t2, 0)
  }

  # Keep top n_secondary zones
  top_idx <- order(-zone_llr)
  top_idx <- top_idx[seq_len(min(n_secondary, sum(zone_llr > 0)))]

  all_pairs <- data.frame(
    center     = zone_centers[top_idx],
    n_regions  = zone_sizes[top_idx],
    cases      = zone_cases[top_idx],
    expected   = zone_expected[top_idx],
    population = zone_pop_vec[top_idx],
    llr        = zone_llr[top_idx],
    stringsAsFactors = FALSE
  )

  result <- list(
    most_likely_cluster = list(
      zone_idx = best_idx,
      region_ids = regions$region_id[best_idx],
      center = zones[[best_z]]$center,
      cases = best_cases,
      expected = expected,
      population = best_pop,
      rr = (best_cases / best_pop) / (C / N),
      llr = obs_llr
    ),
    pvalue = pvalue,
    alpha = alpha,
    nsim = nsim,
    total_cases = C,
    total_population = N,
    regions = regions,
    simulated_llr = sim_llr,
    secondary_clusters = all_pairs,
    model = model
  )

  class(result) <- "circular_scan"
  result
}
