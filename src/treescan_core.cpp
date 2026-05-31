// treescan_core.cpp
// Core C++ functions for the treeSS package
//
// Zones stored in CSR (compressed sparse row) format for cache efficiency.
// Fused multiply + LLR + max in a single pass (no temp matrix allocation).
//
// Monte Carlo simulation uses a single native implementation for all
// n_cores >= 1. Each simulation draws from a deterministic per-simulation
// seed (taken from R's RNG when `seed` is set), so the simulated null
// distribution is identical for any n_cores: n_cores changes only how many
// OpenMP threads share the work, never the result. A thread-local
// std::mt19937 is used because R's RNG is not thread-safe. This replaces
// the earlier design, in which n_cores = 1 used R's rmultinom on a separate
// code path; serial and parallel timings are now directly comparable.
//
// Supports Poisson and Binomial models (model: 0 = Poisson, 1 = Binomial).

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <random>
#include <cstdint>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// --- LLR functions ---

// Poisson LLR for a single (branch, zone) pair
static inline double poisson_llr(double cz, double nz,
                                  double C_branch, double N) {
  double expected = C_branch * nz / N;
  if (cz <= expected || cz <= 0.0 || expected <= 0.0) return 0.0;

  double cz_bar = C_branch - cz;
  double exp_bar = C_branch - expected;

  double llr = cz * std::log(cz / expected);
  if (cz_bar > 0.0 && exp_bar > 0.0) {
    llr += cz_bar * std::log(cz_bar / exp_bar);
  }
  return llr;
}

// Binomial LLR for a single (branch, zone) pair
// cz = cases in zone for branch, nz = population in zone
// C_branch = total cases for branch, N = total population
static inline double binomial_llr(double cz, double nz,
                                   double C_branch, double N) {
  if (nz <= 0.0 || N <= 0.0) return 0.0;

  double nz_bar = N - nz;
  double cz_bar = C_branch - cz;

  // Only consider if rate inside > rate outside
  if (nz_bar <= 0.0) return 0.0;
  if (cz / nz <= cz_bar / nz_bar) return 0.0;
  if (cz <= 0.0) return 0.0;

  // Binomial log-likelihood ratio
  double llr = 0.0;

  // Inside zone
  if (cz > 0.0 && cz < nz) {
    llr += cz * std::log(cz / nz) + (nz - cz) * std::log(1.0 - cz / nz);
  } else if (cz >= nz) {
    llr += cz * std::log(cz / nz);
  }

  // Outside zone
  if (cz_bar > 0.0 && cz_bar < nz_bar) {
    llr += cz_bar * std::log(cz_bar / nz_bar) +
           (nz_bar - cz_bar) * std::log(1.0 - cz_bar / nz_bar);
  } else if (cz_bar >= nz_bar) {
    llr += cz_bar * std::log(cz_bar / nz_bar);
  }

  // Under H0
  if (C_branch > 0.0 && C_branch < N) {
    llr -= C_branch * std::log(C_branch / N) +
           (N - C_branch) * std::log(1.0 - C_branch / N);
  } else if (C_branch >= N) {
    llr -= C_branch * std::log(C_branch / N);
  }

  return (llr > 0.0) ? llr : 0.0;
}

// Dispatch LLR by model
static inline double compute_llr(double cz, double nz,
                                  double C_branch, double N, int model) {
  if (model == 1) {
    return binomial_llr(cz, nz, C_branch, N);
  } else {
    return poisson_llr(cz, nz, C_branch, N);
  }
}

// --- Thread-safe multinomial sampling ---
//
// Conditional binomial method: to sample (X_1, ..., X_k) ~ Multinomial(n, p),
// sample X_1 ~ Bin(n, p_1), then X_2 ~ Bin(n - X_1, p_2 / (1 - p_1)), ...
// This is exact and thread-safe if the generator is thread-local.

static void rmultinom_thread(int n,
                              const double* prob,
                              int k,
                              int* out,
                              std::mt19937& gen) {
  int remaining = n;
  double p_rest = 1.0;

  for (int i = 0; i < k; i++) {
    if (remaining <= 0 || p_rest <= 0.0) {
      out[i] = 0;
      continue;
    }
    if (i == k - 1) {
      out[i] = remaining;
      break;
    }
    double p_i = prob[i];
    double p_cond = p_i / p_rest;
    if (p_cond >= 1.0) p_cond = 1.0;
    if (p_cond <= 0.0) {
      out[i] = 0;
    } else {
      std::binomial_distribution<int> binom(remaining, p_cond);
      int x = binom(gen);
      out[i] = x;
      remaining -= x;
    }
    p_rest -= p_i;
    if (p_rest < 0.0) p_rest = 0.0;
  }
}

// --- Aggregate tree bottom-up ---

// aggregate_up variant that works on flat vector (for thread-local buffers)
static void aggregate_up_flat(double* mat, int n_rows, int n_cols,
                               const int* children_idx,
                               const int* children_ptr,
                               const int* proc_order, int n_order) {
  for (int oi = 0; oi < n_order; oi++) {
    int node = proc_order[oi];
    int ch_start = children_ptr[node];
    int ch_end   = children_ptr[node + 1];
    if (ch_end > ch_start) {
      for (int j = 0; j < n_cols; j++) {
        double s = 0.0;
        for (int c = ch_start; c < ch_end; c++) {
          s += mat[children_idx[c] * n_cols + j];
        }
        mat[node * n_cols + j] = s;
      }
    }
  }
}

// --- Max LLR across all (node, zone) pairs ---

// Flat version for thread-local buffers
static double max_llr_all_pairs_flat(const double* full_cases,
                                      int n_nodes, int n_regions,
                                      const double* Cg, double N,
                                      const int* zone_region_idx,
                                      const int* zone_ptr,
                                      const double* zone_pop, int n_zones,
                                      int model) {
  double max_llr = 0.0;
  for (int g = 0; g < n_nodes; g++) {
    double Cb = Cg[g];
    if (Cb <= 0.0) continue;
    const double* row_g = full_cases + (size_t)g * n_regions;
    for (int z = 0; z < n_zones; z++) {
      double cz = 0.0;
      int zstart = zone_ptr[z];
      int zend   = zone_ptr[z + 1];
      for (int p = zstart; p < zend; p++) {
        cz += row_g[zone_region_idx[p]];
      }
      double llr = compute_llr(cz, zone_pop[z], Cb, N, model);
      if (llr > max_llr) max_llr = llr;
    }
  }
  return max_llr;
}


// ========================================================================
// Monte Carlo for tree-spatial scan statistic
// model: 0 = Poisson, 1 = Binomial
// n_cores: number of OpenMP threads (>=1); the result is independent of it
// ========================================================================
// [[Rcpp::export]]
List mc_treespatial_cpp(NumericMatrix full_cases,
                        NumericVector Cg,
                        double N,
                        IntegerVector zone_region_idx,
                        IntegerVector zone_ptr,
                        NumericVector zone_pop,
                        IntegerVector leaf_tree_idx,
                        IntegerVector children_idx,
                        IntegerVector children_ptr,
                        IntegerVector proc_order,
                        NumericVector prob,
                        int nsim,
                        int model = 0,
                        int n_cores = 1) {

  int n_nodes   = full_cases.nrow();
  int n_regions = full_cases.ncol();
  int n_zones   = zone_pop.size();
  int n_leaves  = leaf_tree_idx.size();

  // Observed data: find best (g, z)
  double obs_max_llr = 0.0;
  int obs_best_g = 0, obs_best_z = 0;

  for (int g = 0; g < n_nodes; g++) {
    double Cb = Cg[g];
    if (Cb <= 0.0) continue;
    for (int z = 0; z < n_zones; z++) {
      double cz = 0.0;
      int zstart = zone_ptr[z];
      int zend   = zone_ptr[z + 1];
      for (int p = zstart; p < zend; p++) {
        cz += full_cases(g, zone_region_idx[p]);
      }
      double llr = compute_llr(cz, zone_pop[z], Cb, N, model);
      if (llr > obs_max_llr) {
        obs_max_llr = llr;
        obs_best_g = g;
        obs_best_z = z;
      }
    }
  }

  NumericVector sim_llr(nsim);

  // ---- Simulation path (single implementation for all n_cores) ----
  // The native sampler over flat arrays is used for every n_cores >= 1:
  // with n_cores = 1 it runs single-threaded, with n_cores > 1 the
  // replicates are split across OpenMP threads. Because each simulation
  // draws from a deterministic per-simulation seed (below), the result is
  // identical for any n_cores given a fixed `seed`; n_cores changes only
  // wall-clock time, never the simulated null distribution. This makes the
  // serial and parallel timings directly comparable.
  //
  // Draw one 32-bit seed per simulation from R's RNG so results are
  // reproducible when `seed` is set in R, independent of n_cores.
  std::vector<uint32_t> seeds(nsim);
  for (int s = 0; s < nsim; s++) {
    // Uniform on [0, 2^32) via R's RNG
    double u = ::unif_rand();
    seeds[s] = (uint32_t)(u * 4294967295.0);
  }

  // Raw C arrays of inputs so we can read them thread-safely
  const double* prob_ptr            = prob.begin();
  const int*    zone_region_idx_ptr = zone_region_idx.begin();
  const int*    zone_ptr_ptr        = zone_ptr.begin();
  const double* zone_pop_ptr        = zone_pop.begin();
  const int*    children_idx_ptr   = children_idx.begin();
  const int*    children_ptr_ptr    = children_ptr.begin();
  const int*    proc_order_ptr      = proc_order.begin();
  const double* Cg_ptr              = Cg.begin();
  const int*    leaf_tree_idx_ptr   = leaf_tree_idx.begin();
  const int     n_order             = proc_order.size();

  double* sim_llr_out = sim_llr.begin();

#ifdef _OPENMP
  omp_set_num_threads(n_cores);
#endif

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
    // Thread-local buffers
    std::vector<double> sim_full((size_t)n_nodes * n_regions, 0.0);
    std::vector<int>    rN(n_regions);
    std::vector<double> sim_Cg(n_nodes);

#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (int s = 0; s < nsim; s++) {
      std::mt19937 gen(seeds[s]);

      std::fill(sim_full.begin(), sim_full.end(), 0.0);

      for (int li = 0; li < n_leaves; li++) {
        int tr = leaf_tree_idx_ptr[li];
        int Cg_leaf = (int)Cg_ptr[tr];
        if (Cg_leaf <= 0) continue;
        rmultinom_thread(Cg_leaf, prob_ptr, n_regions,
                         rN.data(), gen);
        double* row_tr = sim_full.data() + (size_t)tr * n_regions;
        for (int j = 0; j < n_regions; j++) {
          row_tr[j] = (double)rN[j];
        }
      }

      aggregate_up_flat(sim_full.data(), n_nodes, n_regions,
                        children_idx_ptr, children_ptr_ptr,
                        proc_order_ptr, n_order);

      for (int g = 0; g < n_nodes; g++) {
        double rs = 0.0;
        const double* row_g = sim_full.data() + (size_t)g * n_regions;
        for (int j = 0; j < n_regions; j++) rs += row_g[j];
        sim_Cg[g] = rs;
      }

      sim_llr_out[s] = max_llr_all_pairs_flat(sim_full.data(),
                                               n_nodes, n_regions,
                                               sim_Cg.data(), N,
                                               zone_region_idx_ptr,
                                               zone_ptr_ptr, zone_pop_ptr,
                                               n_zones, model);
    }
  }

  return List::create(
    Named("obs_max_llr") = obs_max_llr,
    Named("obs_best_g")  = obs_best_g + 1,
    Named("obs_best_z")  = obs_best_z + 1,
    Named("sim_llr")     = sim_llr
  );
}


// ========================================================================
// Monte Carlo for circular spatial scan statistic
// model: 0 = Poisson, 1 = Binomial
// ========================================================================
// [[Rcpp::export]]
List mc_spatial_cpp(NumericVector cases,
                    double N,
                    double C,
                    IntegerVector zone_region_idx,
                    IntegerVector zone_ptr,
                    NumericVector zone_pop,
                    NumericVector prob,
                    int nsim,
                    int model = 0,
                    int n_cores = 1) {

  int n_regions = cases.size();
  int n_zones   = zone_pop.size();

  double obs_max_llr = 0.0;
  int obs_best_z = 0;

  for (int z = 0; z < n_zones; z++) {
    double cz = 0.0;
    int zstart = zone_ptr[z];
    int zend   = zone_ptr[z + 1];
    for (int p = zstart; p < zend; p++) {
      cz += cases[zone_region_idx[p]];
    }
    double llr = compute_llr(cz, zone_pop[z], C, N, model);
    if (llr > obs_max_llr) {
      obs_max_llr = llr;
      obs_best_z = z;
    }
  }

  NumericVector sim_llr(nsim);

  // ---- Simulation path (single implementation for all n_cores) ----
  // See mc_treespatial_cpp: one native implementation for every
  // n_cores >= 1, identical results for a fixed seed regardless of
  // n_cores, so serial and parallel timings are directly comparable.
  std::vector<uint32_t> seeds(nsim);
  for (int s = 0; s < nsim; s++) {
    double u = ::unif_rand();
    seeds[s] = (uint32_t)(u * 4294967295.0);
  }

  const double* prob_ptr            = prob.begin();
  const int*    zone_region_idx_ptr = zone_region_idx.begin();
  const int*    zone_ptr_ptr        = zone_ptr.begin();
  const double* zone_pop_ptr        = zone_pop.begin();
  double*       sim_llr_out         = sim_llr.begin();

#ifdef _OPENMP
  omp_set_num_threads(n_cores);
#endif

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
    std::vector<int> rN(n_regions);

#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (int s = 0; s < nsim; s++) {
      std::mt19937 gen(seeds[s]);
      rmultinom_thread((int)C, prob_ptr, n_regions, rN.data(), gen);

      double max_llr = 0.0;
      for (int z = 0; z < n_zones; z++) {
        double cz = 0.0;
        int zstart = zone_ptr_ptr[z];
        int zend   = zone_ptr_ptr[z + 1];
        for (int p = zstart; p < zend; p++) {
          cz += (double)rN[zone_region_idx_ptr[p]];
        }
        double llr = compute_llr(cz, zone_pop_ptr[z], C, N, model);
        if (llr > max_llr) max_llr = llr;
      }
      sim_llr_out[s] = max_llr;
    }
  }

  return List::create(
    Named("obs_max_llr") = obs_max_llr,
    Named("obs_best_z")  = obs_best_z + 1,
    Named("sim_llr")     = sim_llr
  );
}


// ========================================================================
// Monte Carlo for tree-based scan statistic
// model: 0 = Poisson, 1 = Binomial
// ========================================================================
// [[Rcpp::export]]
List mc_treescan_cpp(NumericVector node_cases,
                     NumericVector node_pop,
                     double C,
                     double N,
                     IntegerVector leaf_tree_idx,
                     NumericVector leaf_pop,
                     IntegerVector children_idx,
                     IntegerVector children_ptr,
                     IntegerVector proc_order,
                     int nsim,
                     int model = 0,
                     int n_cores = 1) {

  int n_nodes  = node_cases.size();
  int n_leaves = leaf_tree_idx.size();

  NumericVector obs_llr(n_nodes);
  double obs_max_llr = 0.0;
  int obs_best = 0;

  for (int i = 0; i < n_nodes; i++) {
    obs_llr[i] = compute_llr(node_cases[i], node_pop[i], C, N, model);
    if (obs_llr[i] > obs_max_llr) {
      obs_max_llr = obs_llr[i];
      obs_best = i;
    }
  }

  NumericVector sim_llr(nsim);

  NumericVector lprob(n_leaves);
  double total_leaf_pop = 0.0;
  for (int i = 0; i < n_leaves; i++) total_leaf_pop += leaf_pop[i];
  for (int i = 0; i < n_leaves; i++) lprob[i] = leaf_pop[i] / total_leaf_pop;

  // ---- Simulation path (single implementation for all n_cores) ----
  // See mc_treespatial_cpp: one native implementation for every
  // n_cores >= 1, identical results for a fixed seed regardless of
  // n_cores, so serial and parallel timings are directly comparable.
  std::vector<uint32_t> seeds(nsim);
  for (int s = 0; s < nsim; s++) {
    double u = ::unif_rand();
    seeds[s] = (uint32_t)(u * 4294967295.0);
  }

  const double* lprob_ptr         = lprob.begin();
  const int*    leaf_tree_idx_ptr = leaf_tree_idx.begin();
  const int*    children_idx_ptr  = children_idx.begin();
  const int*    children_ptr_ptr   = children_ptr.begin();
  const int*    proc_order_ptr    = proc_order.begin();
  const double* node_pop_ptr      = node_pop.begin();
  double*       sim_llr_out       = sim_llr.begin();
  const int     n_order           = proc_order.size();

#ifdef _OPENMP
  omp_set_num_threads(n_cores);
#endif

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
    std::vector<double> sim_node_cases(n_nodes);
    std::vector<int>    rN(n_leaves);

#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (int s = 0; s < nsim; s++) {
      std::mt19937 gen(seeds[s]);
      std::fill(sim_node_cases.begin(), sim_node_cases.end(), 0.0);

      rmultinom_thread((int)C, lprob_ptr, n_leaves, rN.data(), gen);
      for (int li = 0; li < n_leaves; li++) {
        sim_node_cases[leaf_tree_idx_ptr[li]] = (double)rN[li];
      }

      for (int oi = 0; oi < n_order; oi++) {
        int node = proc_order_ptr[oi];
        int ch_start = children_ptr_ptr[node];
        int ch_end   = children_ptr_ptr[node + 1];
        if (ch_end > ch_start) {
          double s_val = 0.0;
          for (int c = ch_start; c < ch_end; c++) {
            s_val += sim_node_cases[children_idx_ptr[c]];
          }
          sim_node_cases[node] = s_val;
        }
      }

      double max_llr = 0.0;
      for (int i = 0; i < n_nodes; i++) {
        double llr = compute_llr(sim_node_cases[i], node_pop_ptr[i],
                                  C, N, model);
        if (llr > max_llr) max_llr = llr;
      }
      sim_llr_out[s] = max_llr;
    }
  }

  return List::create(
    Named("obs_llr")     = obs_llr,
    Named("obs_max_llr") = obs_max_llr,
    Named("obs_best")    = obs_best + 1,
    Named("sim_llr")     = sim_llr
  );
}

