// src/loglik_math.cpp
//
// R-callable pure-C++ math-scale log-likelihood. Same inputs as the R
// function loglik_math, same semantics (including the `negative` flag
// and mask handling). Intended both as a user-facing drop-in alternative
// to R::loglik_math and as the parity reference for the XPtr hot path.

#include <Rcpp.h>
#include <vector>
#include <string>

#include "loglik_math.h"

// Forward-declare the canonical-vector assembler from math_to_bio.cpp.
Rcpp::NumericVector build_canonical_param_vector_cpp(
    Rcpp::NumericVector param_vector,
    Rcpp::Nullable<Rcpp::NumericVector> mask,
    int p
);

//' Pure-C++ log-likelihood for the xsdm model (math-scale parameters)
//'
//' Computes the log-likelihood directly in C++ without any R callback in
//' the inner loop. Semantically equivalent to the R function loglik_math.
//'
//' @param param_vector Named numeric vector of math-scale parameters.
//'   When `mask` is NULL, must contain every canonical name for the
//'   dimension p implied by `env_dat`. When `mask` is supplied, contains
//'   only the free (non-masked) parameters.
//' @param env_dat 3D numeric array with dimensions (n_loc, ts_length, p).
//'   No missing values allowed.
//' @param occ Integer or logical vector of length n_loc, 0/1 or FALSE/TRUE.
//' @param mask Optional named numeric vector of fixed parameters.
//' @param negative Logical; if TRUE (default) returns the negative
//'   log-likelihood (the value to be minimized).
//' @param num_threads Integer; 0 leaves the RcppParallel default.
//'
//' @returns A scalar double.
// [[Rcpp::export]]
double loglik_math_cpp(
    Rcpp::NumericVector param_vector,
    Rcpp::NumericVector env_dat,
    Rcpp::IntegerVector occ,
    Rcpp::Nullable<Rcpp::NumericVector> mask = R_NilValue,
    bool negative = true,
    int  num_threads = 0
) {
    // Extract 3D dims from the env_dat attribute.
    Rcpp::IntegerVector dims;
    if (!env_dat.hasAttribute("dim")) {
        Rcpp::stop("`env_dat` must be a 3D array.");
    }
    dims = env_dat.attr("dim");
    if (dims.size() != 3) {
        Rcpp::stop("`env_dat` must have 3 dimensions.");
    }
    const int n_loc     = dims[0];
    const int ts_length = dims[1];
    const int p         = dims[2];

    if (occ.size() != n_loc) {
        Rcpp::stop("`occ` length must equal n_loc (dim(env_dat)[1]).");
    }
    if (env_dat.size() != static_cast<R_xlen_t>(n_loc) * ts_length * p) {
        Rcpp::stop("`env_dat` size does not match its dim attribute.");
    }

    // Build the full canonical parameter vector once. Using the same
    // validator that the XPtr closure uses keeps error messages consistent.
    Rcpp::NumericVector full = build_canonical_param_vector_cpp(
        param_vector, mask, p
    );

    // Thread management mirrors log_prob_detect_cpp / loglik_bio_cpp.
    Rcpp::Environment rcppPar = Rcpp::Environment::namespace_env("RcppParallel");
    Rcpp::Function defaultNT  = rcppPar["defaultNumThreads"];
    Rcpp::Function setTO      = rcppPar["setThreadOptions"];
    int old_threads = Rcpp::as<int>(defaultNT());
    if (num_threads > 0) setTO(Rcpp::Named("numThreads") = num_threads);

    std::vector<int> occ_std(occ.begin(), occ.end());

    xsdm::BioParams bp;
    bp.resize(p);
    const double ll = xsdm::loglik_math_eval(
        &full[0],
        p,
        &env_dat[0],
        occ_std.data(),
        n_loc, ts_length,
        bp
    );

    setTO(Rcpp::Named("numThreads") = old_threads);

    return negative ? -ll : ll;
}
