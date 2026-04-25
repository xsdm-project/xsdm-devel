// src/math_to_bio.cpp
//
// R-callable wrapper for xsdm::math_to_bio_apply and the pure-C++ version
// of create_param_vector_masked(). Both are internal (prefixed with `.`) —
// they exist so tests can assert parity with the R reference implementation,
// and so the XPtr closure in loglik_math_xptr.cpp has a single place to
// construct a canonical parameter vector from (free-params + mask).

#include <Rcpp.h>
#include <algorithm>
#include <string>
#include <vector>

#include "math_to_bio.hpp"

// ---------------------------------------------------------------------------
// Canonical parameter-name schema, matching R/make_mask_names.R:
//   mu1..mup, sigltil1..sigltilp, sigrtil1..sigrtilp, ctil, pd, o_par1..o_parq
// ---------------------------------------------------------------------------
static std::vector<std::string> canonical_names_cpp(int p) {
    std::vector<std::string> out;
    out.reserve(xsdm::num_par_cpp(p));
    for (int i = 1; i <= p; ++i) out.push_back("mu" + std::to_string(i));
    for (int i = 1; i <= p; ++i) out.push_back("sigltil" + std::to_string(i));
    for (int i = 1; i <= p; ++i) out.push_back("sigrtil" + std::to_string(i));
    out.push_back("ctil");
    out.push_back("pd");
    const int q = p * (p - 1) / 2;
    for (int i = 1; i <= q; ++i) out.push_back("o_par" + std::to_string(i));
    return out;
}

// Infer p from total number of parameters via num_par(p) = 3p + 2 + p(p-1)/2.
static int infer_p_from_n(int n) {
    if (n < 5) Rcpp::stop("parameter vector too short; need at least 5.");
    const double disc = 9.0 + 8.0 * static_cast<double>(n);
    const double sq   = std::sqrt(disc);
    const double p_d  = (-5.0 + sq) / 2.0;
    const int p = static_cast<int>(std::round(p_d));
    if (p < 1 || xsdm::num_par_cpp(p) != n) {
        Rcpp::stop("Invalid parameter-vector length: %d does not match any num_par(p).", n);
    }
    return p;
}

// Same contract as the R function math_to_bio: expects a named numeric vector
// in canonical order (length == num_par(p)) and returns an R list with
// elements mu, sigltil, sigrtil, ctil, pd, o_mat.

// [[Rcpp::export(.math_to_bio_cpp)]]
Rcpp::List math_to_bio_cpp(Rcpp::NumericVector param_vector) {
    if (!param_vector.hasAttribute("names")) {
        Rcpp::stop("`param_vector` must be a named numeric vector.");
    }
    const int n = param_vector.size();
    const int p = infer_p_from_n(n);

    // Verify canonical ordering of names.
    Rcpp::CharacterVector got_names = param_vector.names();
    const auto expected = canonical_names_cpp(p);
    for (int i = 0; i < n; ++i) {
        if (std::string(got_names[i]) != expected[i]) {
            Rcpp::stop(
                "`param_vector` names do not match canonical order at position %d: "
                "expected '%s', got '%s'.",
                i + 1, expected[i].c_str(), std::string(got_names[i]).c_str()
            );
        }
    }

    xsdm::BioParams bp;
    xsdm::math_to_bio_apply(&param_vector[0], p, bp);

    // Pack o_mat as p x p R matrix (column-major, same layout we wrote).
    Rcpp::NumericMatrix o_mat(p, p);
    std::copy(bp.o_mat.begin(), bp.o_mat.end(), &o_mat[0]);

    return Rcpp::List::create(
        Rcpp::Named("mu")      = Rcpp::NumericVector(bp.mu.begin(), bp.mu.end()),
        Rcpp::Named("sigltil") = Rcpp::NumericVector(bp.sigltil.begin(), bp.sigltil.end()),
        Rcpp::Named("sigrtil") = Rcpp::NumericVector(bp.sigrtil.begin(), bp.sigrtil.end()),
        Rcpp::Named("ctil")    = bp.ctil,
        Rcpp::Named("pd")      = bp.pd,
        Rcpp::Named("o_mat")   = o_mat
    );
}

// Pure C++ version of create_param_vector_masked from R. Merges the optional
// mask with param_vector and returns a fully-populated named numeric vector
// in canonical order. Errors out on unknown names, missing slots, or overlap.
// Used by the XPtr closure in loglik_math_xptr.cpp.

// [[Rcpp::export(.build_canonical_param_vector_cpp)]]
Rcpp::NumericVector build_canonical_param_vector_cpp(
    Rcpp::NumericVector param_vector,
    Rcpp::Nullable<Rcpp::NumericVector> mask,
    int p
) {
    const auto expected = canonical_names_cpp(p);
    const int N = static_cast<int>(expected.size());

    Rcpp::NumericVector out(N, NA_REAL);
    out.attr("names") = Rcpp::wrap(expected);

    auto idx_of = [&expected](const std::string& nm) -> int {
        auto it = std::find(expected.begin(), expected.end(), nm);
        return it == expected.end() ? -1 : static_cast<int>(it - expected.begin());
    };

    if (mask.isNotNull()) {
        Rcpp::NumericVector m(mask);
        Rcpp::CharacterVector m_names = m.names();
        for (int i = 0; i < m.size(); ++i) {
            const std::string nm(m_names[i]);
            const int k = idx_of(nm);
            if (k < 0) Rcpp::stop("Unexpected name in `mask`: '%s'.", nm.c_str());
            out[k] = m[i];
        }
    }

    if (!param_vector.hasAttribute("names")) {
        Rcpp::stop("`param_vector` must be a named numeric vector.");
    }
    Rcpp::CharacterVector pv_names = param_vector.names();
    for (int i = 0; i < param_vector.size(); ++i) {
        const std::string nm(pv_names[i]);
        const int k = idx_of(nm);
        if (k < 0) Rcpp::stop("Unexpected name in `param_vector`: '%s'.", nm.c_str());
        if (!Rcpp::NumericVector::is_na(out[k]) && mask.isNotNull()) {
            Rcpp::NumericVector m(mask);
            Rcpp::CharacterVector m_names = m.names();
            for (int j = 0; j < m.size(); ++j) {
                if (std::string(m_names[j]) == nm) {
                    Rcpp::stop("`param_vector` and `mask` overlap on '%s'.", nm.c_str());
                }
            }
        }
        out[k] = param_vector[i];
    }

    for (int i = 0; i < N; ++i) {
        if (Rcpp::NumericVector::is_na(out[i])) {
            Rcpp::stop("Missing canonical parameter: '%s'.", expected[i].c_str());
        }
    }
    return out;
}
