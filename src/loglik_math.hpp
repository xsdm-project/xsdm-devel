// src/loglik_math.hpp
//
// Pure-C++ math-scale log-likelihood for the xsdm model. Composed of
// math_to_bio_apply + loglik_bio_tile. The R-callable wrapper is in
// loglik_math.cpp; this header exposes the zero-allocation hot-path
// entry used by the XPtr ObjFun closure in loglik_math_xptr.cpp.

#pragma once

#include <vector>

#include "loglik_bio.hpp"
#include "math_to_bio.hpp"

namespace xsdm {

// Zero-allocation hot path. All allocations (bp resize, log_p scratch)
// happen inside loglik_bio_tile which the optimizer calls many times —
// the BioParams itself is reused across calls by the closure.
//
// param_ptr:    length num_par(p), canonical order
// env_dat_ptr:  column-major flat, (n_loc, ts_length, p)
// occ_ptr:      length n_loc, 0 or 1
// bp_scratch:   a BioParams owned by the caller; this function
//               overwrites its fields. Passing it in avoids a vector
//               re-allocation on every optimizer call.
inline double loglik_math_eval(
    const double* param_ptr,
    int           p,
    const double* env_dat_ptr,
    const int*    occ_ptr,
    int           n_loc,
    int           ts_length,
    BioParams&    bp_scratch
) {
    math_to_bio_apply(param_ptr, p, bp_scratch);
    return loglik_bio_tile(env_dat_ptr, occ_ptr, n_loc, ts_length, p, bp_scratch);
}

} // namespace xsdm
