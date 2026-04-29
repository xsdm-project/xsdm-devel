// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

#include "log_prob_detect.h"

#include <algorithm>
#include <cmath>
#include <vector>

using namespace RcppParallel;

// ---------------------------------------------------------------------------
// LPDColumnWorker
//
// RcppParallel Worker that computes the per-(location, time) sum-of-squared
// asymmetric projections.  Replicates the ColumnWorker in like_ltsg.cpp, but
// accesses env_dat_vec directly (avoids building the permuted matrix) using
// the index mapping:
//
//   env_dat_mat[k, j] = env_dat_ptr[l + n_loc*t + n_loc*ts_length*k]
//
// where j = t + ts_length * l  (time t varies fastest, so consecutive j's
// share the same location).
//
// Asymmetric scaling (matches like_ltsg.cpp):
//   usym_i   = inv_sigl[i] * dot_i + inv_drl[i] * max(0, dot_i)
//   col_sum += usym_i^2
// ---------------------------------------------------------------------------
struct LPDColumnWorker : public Worker {
    const double* env_dat_ptr;
    int n_loc;
    int ts_length;
    int p;
    const double* inv_sigl;   // length p: 1/sigltil[i]
    const double* inv_drl;    // length p: 1/sigrtil[i] - 1/sigltil[i]
    const double* o_mat_ptr;  // column-major p x p; o_mat_ptr[k + p*i] = o_mat[k,i] = t(o_mat)[i,k]
    const double* mu;         // length p
    double*       output;     // length ts_length * n_loc

    LPDColumnWorker(
        const double* env_dat_ptr_,
        int n_loc_,
        int ts_length_,
        int p_,
        const double* inv_sigl_,
        const double* inv_drl_,
        const double* o_mat_ptr_,
        const double* mu_,
        double* output_
    )
      : env_dat_ptr(env_dat_ptr_),
        n_loc(n_loc_),
        ts_length(ts_length_),
        p(p_),
        inv_sigl(inv_sigl_),
        inv_drl(inv_drl_),
        o_mat_ptr(o_mat_ptr_),
        mu(mu_),
        output(output_)
    {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t j = begin; j < end; j++) {
            // Column j encodes time t and location l: j = t + ts_length * l
            int l = static_cast<int>(j) / ts_length;
            int t = static_cast<int>(j) % ts_length;

            double col_sum = 0.0;

            for (int i = 0; i < p; i++) {
                double dot = 0.0;

                for (int k = 0; k < p; k++) {
                    double env_val = env_dat_ptr[
                        l +
                        static_cast<long long>(n_loc) * t +
                        static_cast<long long>(n_loc) * ts_length * k
                    ];

                    // t(o_mat)[i, k] = o_mat[k, i] = o_mat_ptr[k + p*i] (col-major)
                    double o_val = o_mat_ptr[k + p * i];

                    dot += o_val * (env_val - mu[k]);
                }

                double usym = inv_sigl[i] * dot + inv_drl[i] * std::max(0.0, dot);
                col_sum += usym * usym;
            }

            output[j] = col_sum;
        }
    }
};

// ---------------------------------------------------------------------------
// xsdm::log_prob_detect_tile  — implementation
// ---------------------------------------------------------------------------
namespace xsdm {

std::vector<double> log_prob_detect_tile(
    const double* env_dat_ptr,
    int n_loc,
    int ts_length,
    int p,
    const double* mu_ptr,
    const double* sigltil_ptr,
    const double* sigrtil_ptr,
    const double* o_mat_ptr,
    double ctil,
    double pd,
    bool return_prob
) {
    // ---- 1. Build per-direction inverse scaling vectors ----
    // Matches like_neg_ltsgr_cpp.R:
    //   dl_inv  = diag(1 / sigltil)
    //   drl_inv = diag(1 / sigrtil) - diag(1 / sigltil)
    // Inf entries in sigltil/sigrtil are handled by IEEE 754 (1/Inf == 0).
    std::vector<double> inv_sigl(p);
    std::vector<double> inv_drl(p);
    for (int i = 0; i < p; i++) {
        double isl = std::isinf(sigltil_ptr[i]) ? 0.0 : 1.0 / sigltil_ptr[i];
        double isr = std::isinf(sigrtil_ptr[i]) ? 0.0 : 1.0 / sigrtil_ptr[i];
        inv_sigl[i] = isl;
        inv_drl[i]  = isr - isl;
    }

    // ---- 2. Parallel LTSG kernel (RcppParallel ColumnWorker) ----
    // col_sums[j] = col_sum for column j = t + ts_length * l,
    // i.e. col_sums[l * ts_length + t].
    int total = ts_length * n_loc;
    std::vector<double> col_sums(total, 0.0);

    LPDColumnWorker worker(
        env_dat_ptr, n_loc, ts_length, p,
        inv_sigl.data(), inv_drl.data(), o_mat_ptr,
        mu_ptr, col_sums.data()
    );
    parallelFor(0, static_cast<std::size_t>(total), worker);

    // ---- 3. Aggregate over time using xtensor ----
    // h[l] = sum_t col_sums[l*ts_length + t] / (2*ts_length).
    // xt::adapt with row-major shape {n_loc, ts_length}:
    //   sums_view(l, t) = col_sums[l * ts_length + t]  ✓
    std::vector<std::size_t> shape = {
        static_cast<std::size_t>(n_loc),
        static_cast<std::size_t>(ts_length)
    };
    auto sums_view = xt::adapt(
        col_sums.data(),
        static_cast<std::size_t>(total),
        xt::no_ownership(),
        shape
    );
    auto h_xt = xt::sum(sums_view, {1}) / (2.0 * ts_length);

    // ---- 4. Compute log detection probability ----
    // log_p[l] = log(pd) - log1pexp(ctil + h[l])
    std::vector<double> result(n_loc);
    double log_pd = std::log(pd);
    for (int l = 0; l < n_loc; l++) {
        double h_l = h_xt(l);
        double val = log_pd - log1pexp(ctil + h_l);
        result[l] = return_prob ? std::exp(val) : val;
    }

    return result;
}

} // namespace xsdm

// ---------------------------------------------------------------------------
// log_prob_detect_cpp — Rcpp::export wrapper
// ---------------------------------------------------------------------------

//' Compute log detection probabilities from a flat environmental data vector
//'
//' C++ implementation of \code{log_prob_detect()} that accepts environmental
//' data as a flat numeric vector with explicit dimension metadata.  This
//' signature is designed for block-by-block raster evaluation where each
//' block is passed as a contiguous vector rather than a 3-D R array.
//'
//' Collapses the R call chain
//' \code{like_neg_ltsgr_cpp() -> like_ltsg()} into a single xtensor-accelerated
//' C++ function.
//'
//' @param env_dat_vec Numeric vector.  Column-major flat representation of a
//'   3-D array with logical dimensions \code{c(n_loc, ts_length, p)}: variable
//'   \code{k} (1-indexed) occupies positions
//'   \code{(k-1)*n_loc*ts_length + 1} to \code{k*n_loc*ts_length},
//'   and within that block pixels (locations) vary fastest.
//' @param env_dat_dims Integer vector of length 3: \code{c(n_loc, ts_length, p)}.
//' @param mu Numeric vector of length \code{p}.  Optimal environmental values.
//' @param sigltil Numeric vector of length \code{p}.  Positive; \code{Inf}
//'   entries are allowed (treated as zero inverse-scale).
//' @param sigrtil Numeric vector of length \code{p}.  Positive; \code{Inf}
//'   entries are allowed.
//' @param o_mat Numeric matrix, \code{p x p} orthogonal.
//' @param ctil Scalar.  Center of the detection-link function.
//' @param pd Scalar in \code{(0, 1]}.  Maximum probability of detection.
//' @param return_prob Logical.  If \code{TRUE}, return probabilities; if
//'   \code{FALSE} (default) return log-probabilities.
//' @param num_threads Integer.  Number of parallel threads.  \code{0}
//'   (default) uses \code{RcppParallel::defaultNumThreads()}.
//' @return Numeric vector of length \code{n_loc}.
// [[Rcpp::export]]
Rcpp::NumericVector log_prob_detect_cpp(
    Rcpp::NumericVector env_dat_vec,
    Rcpp::IntegerVector env_dat_dims,
    Rcpp::NumericVector mu,
    Rcpp::NumericVector sigltil,
    Rcpp::NumericVector sigrtil,
    Rcpp::NumericMatrix o_mat,
    double ctil,
    double pd,
    bool return_prob = false,
    int  num_threads = 0
) {
    // ---- Input validation ----
    if (env_dat_dims.size() != 3) {
        Rcpp::stop("`env_dat_dims` must have length 3 (n_loc, ts_length, p).");
    }
    int n_loc     = env_dat_dims[0];
    int ts_length = env_dat_dims[1];
    int p         = env_dat_dims[2];

    if (mu.size() != p) {
        Rcpp::stop("`mu` must have length p (= env_dat_dims[3]).");
    }
    if (sigltil.size() != p) {
        Rcpp::stop("`sigltil` must have length p.");
    }
    if (sigrtil.size() != p) {
        Rcpp::stop("`sigrtil` must have length p.");
    }
    if (o_mat.nrow() != p || o_mat.ncol() != p) {
        Rcpp::stop("`o_mat` must be a p x p matrix.");
    }
    if (env_dat_vec.size() != static_cast<R_xlen_t>(n_loc) * ts_length * p) {
        Rcpp::stop("`env_dat_vec` length must equal prod(env_dat_dims).");
    }

    // ---- Thread management: save -> set -> restore ----
    Rcpp::Environment rcppPar  = Rcpp::Environment::namespace_env("RcppParallel");
    Rcpp::Function defaultNT   = rcppPar["defaultNumThreads"];
    Rcpp::Function setTO       = rcppPar["setThreadOptions"];
    int old_threads = Rcpp::as<int>(defaultNT());

    if (num_threads > 0) {
        setTO(Rcpp::Named("numThreads") = num_threads);
    }

    // ---- Delegate to the xtensor-backed tile function ----
    std::vector<double> res = xsdm::log_prob_detect_tile(
        REAL(env_dat_vec),
        n_loc, ts_length, p,
        REAL(mu),
        REAL(sigltil),
        REAL(sigrtil),
        REAL(o_mat),
        ctil, pd, return_prob
    );

    // ---- Restore thread count ----
    setTO(Rcpp::Named("numThreads") = old_threads);

    return Rcpp::NumericVector(res.begin(), res.end());
}
