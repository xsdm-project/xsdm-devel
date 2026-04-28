// src/loglik_bio.cpp
//
// R-callable wrapper around xsdm::loglik_bio_tile. Semantically identical
// to the R function loglik_bio with sum_log_p = TRUE, return_prob = FALSE.
// Other flag combinations (sum_log_p = FALSE, return_prob = TRUE) are not
// needed by the optimizer and intentionally not supported here to keep the
// hot-path contract narrow.

#include <Rcpp.h>
#include <algorithm>
#include <vector>

#include "loglik_bio.h"

//' Pure-C++ log-likelihood for the xsdm model (biological-scale parameters)
//'
//' Computes the log-likelihood directly in C++ without any R callback.
//' Equivalent to loglik_bio(..., sum_log_p = TRUE, return_prob = FALSE).
//'
//' @param env_dat_vec Flat numeric vector containing env_dat in column-major
//'   order (as produced by as.vector(env_dat)).
//' @param env_dat_dims Integer vector of length 3: c(n_loc, ts_length, p).
//' @param occ Integer vector of length n_loc, 0 or 1.
//' @param mu Numeric vector, length p.
//' @param sigltil Positive numeric vector, length p.
//' @param sigrtil Positive numeric vector, length p.
//' @param o_mat A p x p orthogonal matrix (column-major).
//' @param ctil Scalar.
//' @param pd Scalar in (0, 1].
//' @param num_threads Number of threads for the inner xtensor kernel
//'   (0 = RcppParallel default).
//' @return Scalar log-likelihood.
// [[Rcpp::export]]
double loglik_bio_cpp(
    Rcpp::NumericVector env_dat_vec,
    Rcpp::IntegerVector env_dat_dims,
    Rcpp::IntegerVector occ,
    Rcpp::NumericVector mu,
    Rcpp::NumericVector sigltil,
    Rcpp::NumericVector sigrtil,
    Rcpp::NumericMatrix o_mat,
    double              ctil,
    double              pd,
    int                 num_threads = 0
) {
    if (env_dat_dims.size() != 3) Rcpp::stop("env_dat_dims must have length 3.");
    const int n_loc     = env_dat_dims[0];
    const int ts_length = env_dat_dims[1];
    const int p         = env_dat_dims[2];

    if (static_cast<std::size_t>(env_dat_vec.size())
        != static_cast<std::size_t>(n_loc) * ts_length * p) {
        Rcpp::stop("env_dat_vec length does not match product of env_dat_dims.");
    }
    if (occ.size() != n_loc) Rcpp::stop("occ must have length n_loc.");
    if (mu.size() != p || sigltil.size() != p || sigrtil.size() != p) {
        Rcpp::stop("mu / sigltil / sigrtil must have length p.");
    }
    if (o_mat.nrow() != p || o_mat.ncol() != p) {
        Rcpp::stop("o_mat must be p x p.");
    }

    // Thread management: honour caller-requested num_threads via the
    // RcppParallel R namespace, restoring the previous value on return.
    // Matches the pattern in log_prob_detect_cpp.
    Rcpp::Environment rcppPar = Rcpp::Environment::namespace_env("RcppParallel");
    Rcpp::Function    defaultNT = rcppPar["defaultNumThreads"];
    Rcpp::Function    setTO     = rcppPar["setThreadOptions"];
    int old_threads = Rcpp::as<int>(defaultNT());
    if (num_threads > 0) {
        setTO(Rcpp::Named("numThreads") = num_threads);
    }

    xsdm::BioParams bp;
    bp.resize(p);
    std::copy(mu.begin(),      mu.end(),      bp.mu.begin());
    std::copy(sigltil.begin(), sigltil.end(), bp.sigltil.begin());
    std::copy(sigrtil.begin(), sigrtil.end(), bp.sigrtil.begin());
    std::copy(&o_mat[0],       &o_mat[0] + p * p, bp.o_mat.begin());
    bp.ctil = ctil;
    bp.pd   = pd;

    std::vector<int> occ_std(occ.begin(), occ.end());

    const double res = xsdm::loglik_bio_tile(
        &env_dat_vec[0],
        occ_std.data(),
        n_loc, ts_length, p,
        bp
    );

    // Restore thread count before returning.
    setTO(Rcpp::Named("numThreads") = old_threads);
    return res;
}
