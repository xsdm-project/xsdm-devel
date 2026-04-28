// src/expm_skew.cpp
//
// Thin R-callable wrapper around xsdm::build_orthogonal_matrix_cpp.
// Exists primarily as a parity target for tests/testthat/test-expm_skew_cpp.R.
// Not part of the user-facing API — kept internal (no roxygen) and only used
// to prove that the C++ expm matches R's expm::expm on the skew-symmetric
// matrices built by build_orthogonal_matrix().

#include <Rcpp.h>
#include "expm_skew.h"

// [[Rcpp::export(.build_orthogonal_matrix_cpp)]]
Rcpp::NumericMatrix build_orthogonal_matrix_cpp(Rcpp::NumericVector entries) {
    const int q = entries.size();

    // Solve k from k*(k-1)/2 = q.
    int p = 1;
    if (q > 0) {
        const double k_real = 0.5 * (1.0 + std::sqrt(1.0 + 8.0 * q));
        p = static_cast<int>(std::round(k_real));
        if (p < 2 || p * (p - 1) / 2 != q) {
            Rcpp::stop(
                "Invalid 'entries' length: must be a triangular number "
                "n = k*(k-1)/2 for some integer k >= 2. Got n = %d.", q
            );
        }
    }

    Rcpp::NumericMatrix out(p, p);
    xsdm::build_orthogonal_matrix_cpp(
        q == 0 ? nullptr : &entries[0],
        p,
        &out[0]    // column-major
    );
    return out;
}
