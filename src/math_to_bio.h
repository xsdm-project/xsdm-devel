// src/math_to_bio.h
//
// Pure-C++ port of R/math_to_bio.R + (mask + free params) -> canonical vector
// assembly from R/create_param_vector_masked.R.
//
// The R layer names its parameters; the C++ core operates on raw double*
// in canonical order:
//
//   index  0 .. p-1              -> mu[1..p]       (identity)
//   index  p .. 2p-1             -> sigltil        (exp)
//   index  2p .. 3p-1            -> sigrtil        (exp)
//   index  3p                    -> ctil           (identity)
//   index  3p + 1                -> pd             (expit)
//   index  3p + 2 .. 3p + 1 + q  -> o_par          (build_orthogonal_matrix)
//
// where q = p * (p - 1) / 2 and total length = num_par(p) = 3p + 2 + q.

#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

#include "expm_skew.h"

namespace xsdm {

// ---------------------------------------------------------------------------
// Plain-old-data bundle of biological-scale parameters. The XPtr closure
// reuses one of these across evaluations to avoid allocation churn.
// ---------------------------------------------------------------------------
struct BioParams {
    std::vector<double> mu;         // length p
    std::vector<double> sigltil;    // length p
    std::vector<double> sigrtil;    // length p
    double              ctil = 0.0;
    double              pd   = 0.0;
    std::vector<double> o_mat;      // p * p, column-major

    void resize(int p) {
        const std::size_t q = static_cast<std::size_t>(p);
        mu.resize(q);
        sigltil.resize(q);
        sigrtil.resize(q);
        o_mat.resize(q * q);
    }
};

// expit(x) = 1 / (1 + exp(-x)), numerically stable for large |x|.
inline double expit_cpp(double x) {
    if (x >= 0.0) {
        const double z = std::exp(-x);
        return 1.0 / (1.0 + z);
    } else {
        const double z = std::exp(x);
        return z / (1.0 + z);
    }
}

// num_par(p) = 3p + 2 + p(p-1)/2
inline int num_par_cpp(int p) {
    return 3 * p + 2 + (p * (p - 1)) / 2;
}

// ---------------------------------------------------------------------------
// math_to_bio_apply
//
// Writes the biological-scale parameters corresponding to the canonical-order
// math-scale vector param_ptr (length = num_par(p)) into `out`.
// No allocations are performed if `out` already has the right shape.
// ---------------------------------------------------------------------------
inline void math_to_bio_apply(const double* param_ptr, int p, BioParams& out) {
    out.resize(p);
    for (int i = 0; i < p; ++i) {
        out.mu[i]      = param_ptr[i];
        out.sigltil[i] = std::exp(param_ptr[p + i]);
        out.sigrtil[i] = std::exp(param_ptr[2 * p + i]);
    }
    out.ctil = param_ptr[3 * p];
    out.pd   = expit_cpp(param_ptr[3 * p + 1]);

    const int q = p * (p - 1) / 2;
    const double* opar = (q > 0) ? (param_ptr + 3 * p + 2) : nullptr;
    xsdm::build_orthogonal_matrix_cpp(opar, p, out.o_mat.data());
}

} // namespace xsdm
