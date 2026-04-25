#pragma once

#include <cmath>
#include <vector>

// Vendored xtensor-stack headers
#include <xtensor/xarray.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xreducer.hpp>

namespace xsdm {

// ---------------------------------------------------------------------------
// Numerically stable log(1 + exp(x))
// Matches R/log1pexp.R defaults (c0 = -37, c1 = 18, c2 = 33.3)
// Reference: Mächler (2012) "Accurately Computing log(1 - exp(-|a|))"
// ---------------------------------------------------------------------------
inline double log1pexp(double x) {
    if (x <= -37.0) return std::exp(x);
    if (x <= 18.0)  return std::log1p(std::exp(x));
    if (x <= 33.3)  return x + std::exp(-x);
    return x;
}

// ---------------------------------------------------------------------------
// log_prob_detect_tile
//
// Core computation: collapses the R call chain
//   like_neg_ltsgr_cpp() -> like_ltsg()
// into a single C++ function that operates on raw pointers.
//
// Parameters
// ----------
// env_dat_ptr  : column-major flat array, logical dims (n_loc x ts_length x p).
//                env_dat_ptr[l + n_loc*t + n_loc*ts_length*k] = env_dat[l,t,k].
// n_loc        : number of locations
// ts_length    : time-series length
// p            : number of environmental variables
// mu_ptr       : optimal environmental values, length p
// sigltil_ptr  : left-width parameters (positive), length p
// sigrtil_ptr  : right-width parameters (positive), length p
// o_mat_ptr    : column-major p x p orthogonal matrix
// ctil         : detection-link center scalar
// pd           : maximum detection probability in (0, 1]
// return_prob  : if true return probabilities, else return log-probabilities
//
// Returns
// -------
// std::vector<double> of length n_loc
// ---------------------------------------------------------------------------
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
);

} // namespace xsdm
