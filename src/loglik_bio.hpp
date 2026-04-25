// src/loglik_bio.hpp
//
// Pure-C++ log-likelihood for the xsdm model, biological-scale parameters.
// Mirrors R/loglik_bio.R with sum_log_p = TRUE, return_prob = FALSE — the
// only settings that the optimizer cares about. Under the hood it reuses
// xsdm::log_prob_detect_tile (xtensor-backed) to compute the per-location
// log-probability vector, then aggregates
//
//     sum_l occ[l] * log_p[l] + (1 - occ[l]) * log1mexp(-log_p[l])
//
// where log1mexp is the numerically-stable log(1 - exp(-a)) from
// Maechler (2012), identical to R/log1mexp.R.

#pragma once

#include <cmath>
#include <cstdint>
#include <vector>

#include "log_prob_detect.hpp"
#include "math_to_bio.hpp"

namespace xsdm {

// log(1 - exp(-a)) with a >= 0. Uses the two-formula split at log(2)
// recommended by Maechler. Matches R/log1mexp.R.
inline double log1mexp_cpp(double a) {
    // Inputs from loglik_bio_tile below are always -log_p with log_p <= 0,
    // so a >= 0. Preserve R's NaN behaviour for a < 0.
    if (a < 0.0) return std::nan("");
    constexpr double kCutoff = 0.6931471805599453;  // log(2)
    if (a <= kCutoff) {
        return std::log(-std::expm1(-a));
    } else {
        return std::log1p(-std::exp(-a));
    }
}

// ---------------------------------------------------------------------------
// loglik_bio_tile
//
// env_dat_ptr layout: column-major flat, logical dims (n_loc, ts_length, p)
//   env_dat_ptr[l + n_loc * t + n_loc * ts_length * k] = env_dat[l, t, k]
// occ_ptr:    length n_loc, 0 or 1
// bp:         biological-scale parameter bundle
// num_threads:   forwarded to log_prob_detect_tile via parallelFor scheduling
// ---------------------------------------------------------------------------
inline double loglik_bio_tile(
    const double* env_dat_ptr,
    const int*    occ_ptr,
    int           n_loc,
    int           ts_length,
    int           p,
    const BioParams& bp
) {
    // per-location log_p of detection
    std::vector<double> log_p = xsdm::log_prob_detect_tile(
        env_dat_ptr,
        n_loc,
        ts_length,
        p,
        bp.mu.data(),
        bp.sigltil.data(),
        bp.sigrtil.data(),
        bp.o_mat.data(),
        bp.ctil,
        bp.pd,
        /* return_prob = */ false
    );

    double acc = 0.0;
    for (int l = 0; l < n_loc; ++l) {
        const double lp = log_p[l];
        if (occ_ptr[l]) {
            acc += lp;
        } else {
            // occ = 0 => (1 - occ) * log1mexp(-lp). lp <= 0 so -lp >= 0.
            acc += log1mexp_cpp(-lp);
        }
    }
    return acc;
}

} // namespace xsdm
