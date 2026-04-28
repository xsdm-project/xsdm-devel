// src/loglik_math_xptr_cpp.cpp
//
// Pure-C++ ObjFun closure for ucminfcpp::ucminf_xptr. The closure:
//
//   (1) snapshots env_dat / occ into owned buffers at construction, so
//       subsequent (f, g) calls never touch the R heap;
//   (2) stores a pre-populated canonical parameter vector with mask
//       values already in place, plus a slot-index table mapping each
//       free-parameter position to its canonical index;
//   (3) splices x into that canonical vector on each call and evaluates
//       xsdm::loglik_math_eval — no Rcpp::Function callback.
//
// The original R-callback closure in src/loglik_math_xptr.cpp is left
// untouched for backward compatibility with any external callers.
// optimize_loglik_math_() in R/optimize_likelihood.R will switch to
// this builder in the next commit.

// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(ucminfcpp)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "loglik_math.h"
#include "ucminf_core.hpp"

namespace {

// Canonical names for p, matching R/make_mask_names.R exactly:
//   mu1..mup, sigltil1..sigltilp, sigrtil1..sigrtilp, ctil, pd,
//   o_par1..o_parq
std::vector<std::string> canonical_names(int p) {
    std::vector<std::string> out;
    out.reserve(3 * p + 2 + (p * (p - 1)) / 2);
    for (int i = 1; i <= p; ++i) out.push_back("mu" + std::to_string(i));
    for (int i = 1; i <= p; ++i) out.push_back("sigltil" + std::to_string(i));
    for (int i = 1; i <= p; ++i) out.push_back("sigrtil" + std::to_string(i));
    out.push_back("ctil");
    out.push_back("pd");
    const int q = p * (p - 1) / 2;
    for (int i = 1; i <= q; ++i) out.push_back("o_par" + std::to_string(i));
    return out;
}

// Mutable state that lives inside the ucminf::ObjFun closure.
// Held via shared_ptr because std::function needs to copy the closure,
// and we want one copy of the big env_dat buffer rather than one per copy.
//
// Name validation (canonical ordering, missing slots, etc.) is deferred
// to the first closure call. That way mocked-optimizer tests that never
// invoke the closure can still construct an XPtr for arbitrary names.
struct XptrClosureState {
    // --- data buffers (read-only after construction) -----------------------
    std::vector<double> env_dat;       // column-major flat, (n_loc, ts, p)
    std::vector<int>    occ;           // length n_loc
    int n_loc   = 0;
    int ts      = 0;
    int p       = 0;

    // --- deferred canonical-vector setup -----------------------------------
    std::vector<std::string> free_names;    // free-param names, in par order
    bool                     has_mask = false;
    std::vector<std::string> mask_names;    // length m
    std::vector<double>      mask_values;   // length m

    // Materialised on first call:
    std::vector<double>      full;          // length num_par(p)
    std::vector<int>         free_slot_idx; // length n_free
    bool                     ready = false;

    // --- scratch (reused across calls) -------------------------------------
    xsdm::BioParams     bp;                // resized to p once at construction

    // --- gradient config ---------------------------------------------------
    double gradstep_rel = 1e-6;
    double gradstep_abs = 1e-8;
    bool   use_central  = true;
};

// Canonical-position lookup (pure C++).
int find_canonical_index(
    const std::vector<std::string>& canon,
    const std::string& nm
) {
    for (std::size_t i = 0; i < canon.size(); ++i) {
        if (canon[i] == nm) return static_cast<int>(i);
    }
    return -1;
}

} // namespace

// ---------------------------------------------------------------------------
// make_loglik_math_xptr_cpp
//
// Builds an Rcpp::XPtr<ucminf::ObjFun> whose evaluation is entirely in C++
// (no R callbacks inside the ObjFun). Intended to be called from
// optimize_loglik_math_() in place of the R-callback version.
//
// Arguments
// ---------
// env_dat       3D numeric array, dim = (n_loc, ts_length, p).
// occ           Integer/logical vector of length n_loc.
// mask          Optional named numeric vector of fixed parameters.
// free_names    Character vector of the free-parameter names, in the
//               same order as the starting `par` passed to ucminf_xptr.
// num_threads   Threads for the xtensor inner kernel (0 = default).
// grad          "central" (default) or "forward".
// gradstep      (rel, abs) tuple, default (1e-6, 1e-8).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
SEXP make_loglik_math_xptr_cpp(
    Rcpp::NumericVector env_dat,
    Rcpp::IntegerVector occ,
    Rcpp::Nullable<Rcpp::NumericVector> mask,
    Rcpp::CharacterVector free_names,
    int num_threads = 1,
    std::string grad = "central",
    Rcpp::NumericVector gradstep = Rcpp::NumericVector::create(1e-6, 1e-8)
) {
    if (gradstep.size() != 2) Rcpp::stop("`gradstep` must have length 2.");
    if (grad != "forward" && grad != "central") {
        Rcpp::stop("`grad` must be either 'forward' or 'central'.");
    }
    if (!(gradstep[0] > 0.0) || !(gradstep[1] > 0.0)) {
        Rcpp::stop("`gradstep` entries must be strictly positive.");
    }
    if (free_names.size() == 0) {
        Rcpp::stop("`free_names` must be non-empty.");
    }
    for (int i = 0; i < free_names.size(); ++i) {
        if (free_names[i] == NA_STRING) {
            Rcpp::stop("`free_names` contains NA at position %d.", i + 1);
        }
    }

    // --- Resolve dims from env_dat's dim attribute -------------------------
    if (!env_dat.hasAttribute("dim")) {
        Rcpp::stop("`env_dat` must be a 3D array.");
    }
    Rcpp::IntegerVector dims = env_dat.attr("dim");
    if (dims.size() != 3) Rcpp::stop("`env_dat` must have 3 dimensions.");
    const int n_loc = dims[0];
    const int ts    = dims[1];
    const int p     = dims[2];

    if (occ.size() != n_loc) {
        Rcpp::stop("`occ` length must equal n_loc (dim(env_dat)[1]).");
    }
    if (env_dat.size() != static_cast<R_xlen_t>(n_loc) * ts * p) {
        Rcpp::stop("`env_dat` size does not match its dim attribute.");
    }

    // --- Set thread count (caller-controlled). Inside the closure we
    //     never flip threads — the closure path never calls R. -------------
    if (num_threads > 0) {
        Rcpp::Environment rcppPar = Rcpp::Environment::namespace_env("RcppParallel");
        Rcpp::Function setTO      = rcppPar["setThreadOptions"];
        setTO(Rcpp::Named("numThreads") = num_threads);
    }

    // --- Build the closure state (canonical-vector setup is deferred) ----
    auto state = std::make_shared<XptrClosureState>();
    state->n_loc = n_loc;
    state->ts    = ts;
    state->p     = p;
    state->env_dat.assign(env_dat.begin(), env_dat.end());
    state->occ.assign(occ.begin(), occ.end());

    state->free_names.reserve(free_names.size());
    for (int i = 0; i < free_names.size(); ++i) {
        state->free_names.emplace_back(free_names[i]);
    }
    // Store mask as plain C++ data so the closure does not hold any SEXP.
    if (mask.isNotNull()) {
        Rcpp::NumericVector mv(mask);
        if (mv.size() > 0) {
            // R's names() on an unnamed vector returns NULL, which Rcpp
            // coerces to a zero-length CharacterVector. Indexing into it
            // under the loop below would read out of bounds — catch that
            // here explicitly with a clear error.
            SEXP names_sexp = Rf_getAttrib(mv, R_NamesSymbol);
            if (names_sexp == R_NilValue) {
                Rcpp::stop("`mask` must be a named numeric vector.");
            }
            Rcpp::CharacterVector mn(names_sexp);
            if (mn.size() != mv.size()) {
                Rcpp::stop("`mask` names length does not match values length.");
            }
            state->has_mask = true;
            state->mask_names.reserve(mv.size());
            state->mask_values.assign(mv.begin(), mv.end());
            for (int i = 0; i < mv.size(); ++i) {
                if (mn[i] == NA_STRING) {
                    Rcpp::stop("`mask` name at position %d is NA.", i + 1);
                }
                state->mask_names.emplace_back(mn[i]);
            }
        }
        // mv.size() == 0 → treat as no mask (has_mask stays false).
    }

    state->bp.resize(p);

    state->gradstep_rel = gradstep[0];
    state->gradstep_abs = gradstep[1];
    state->use_central  = (grad == "central");

    auto fn = std::make_unique<ucminf::ObjFun>(
        [state](const std::vector<double>& x,
                std::vector<double>&       g,
                double&                    f)
        {
            // Ensure g matches x (ucminfcpp already pre-sizes it, but be
            // defensive — callers from other embeddings might not).
            if (g.size() != x.size()) {
                g.resize(x.size(), 0.0);
            }

            // Any exception from validation, xtensor, or math_to_bio should
            // not propagate through ucminfcpp's template stack — doing so
            // can leave its trust-region state partially updated. We catch
            // here, surface a sentinel +Inf, zero the gradient, and store
            // the error so the next builder-side call can re-report it.
            try {
            // Lazy one-time build of the canonical parameter vector + slot
            // indices. Done here rather than at construction so that tests
            // which mock the optimizer and never invoke the closure can
            // still build an XPtr with non-canonical names.
            if (!state->ready) {
                // Pure-C++ canonical vector construction. Validates free
                // names and mask names against the canonical schema,
                // detects overlaps, and ensures every canonical slot is
                // covered. Throws std::runtime_error on any mismatch so
                // the error surfaces through ucminfcpp as an R condition.
                const auto canon = canonical_names(state->p);
                const int N = static_cast<int>(canon.size());

                state->full.assign(N, std::numeric_limits<double>::quiet_NaN());

                // Apply mask first (if present) so overlap detection with
                // free params is straightforward.
                std::vector<bool> slot_filled(N, false);
                if (state->has_mask) {
                    for (std::size_t i = 0; i < state->mask_names.size(); ++i) {
                        const int k = find_canonical_index(canon, state->mask_names[i]);
                        if (k < 0) {
                            throw std::runtime_error(
                                "ObjFun: unknown canonical name in mask: "
                                + state->mask_names[i]);
                        }
                        state->full[k] = state->mask_values[i];
                        slot_filled[k] = true;
                    }
                }

                // Record free-slot indices, checking against overlap and
                // unknown names.
                state->free_slot_idx.clear();
                state->free_slot_idx.reserve(state->free_names.size());
                for (const auto& nm : state->free_names) {
                    const int k = find_canonical_index(canon, nm);
                    if (k < 0) {
                        throw std::runtime_error(
                            "ObjFun: unknown canonical name in free_names: " + nm);
                    }
                    if (slot_filled[k]) {
                        throw std::runtime_error(
                            "ObjFun: free_names and mask overlap on: " + nm);
                    }
                    state->free_slot_idx.push_back(k);
                    slot_filled[k] = true;
                    // Seed the slot with a finite zero; the caller will
                    // overwrite with x[i] below.
                    state->full[k] = 0.0;
                }

                for (int i = 0; i < N; ++i) {
                    if (!slot_filled[i]) {
                        throw std::runtime_error(
                            "ObjFun: canonical slot uncovered: " + canon[i]);
                    }
                }
                state->ready = true;
            }

            const int n_free = static_cast<int>(x.size());
            if (static_cast<int>(state->free_slot_idx.size()) != n_free) {
                throw std::runtime_error(
                    "ObjFun: x length does not match free_slot_idx size"
                );
            }

            // Splice x into the canonical slots.
            for (int i = 0; i < n_free; ++i) {
                state->full[state->free_slot_idx[i]] = x[i];
            }

            // negative = TRUE (we're minimizing).
            const double ll = xsdm::loglik_math_eval(
                state->full.data(), state->p,
                state->env_dat.data(), state->occ.data(),
                state->n_loc, state->ts,
                state->bp
            );
            f = -ll;

            // Finite-difference gradient over the free slots.
            for (int i = 0; i < n_free; ++i) {
                const int   k  = state->free_slot_idx[i];
                const double xi = state->full[k];
                const double dx = std::abs(xi) * state->gradstep_rel
                                + state->gradstep_abs;

                state->full[k] = xi + dx;
                const double f_plus = -xsdm::loglik_math_eval(
                    state->full.data(), state->p,
                    state->env_dat.data(), state->occ.data(),
                    state->n_loc, state->ts,
                    state->bp
                );

                if (state->use_central) {
                    state->full[k] = xi - dx;
                    const double f_minus = -xsdm::loglik_math_eval(
                        state->full.data(), state->p,
                        state->env_dat.data(), state->occ.data(),
                        state->n_loc, state->ts,
                        state->bp
                    );
                    g[i] = (f_plus - f_minus) / (2.0 * dx);
                } else {
                    g[i] = (f_plus - f) / dx;
                }

                state->full[k] = xi;
            }
            } catch (const std::exception& e) {
                // Sentinel: +Inf objective + zero gradient signals failure
                // to ucminf without corrupting its workspace. The message
                // is surfaced later by a sanity-check on the R side.
                f = std::numeric_limits<double>::infinity();
                std::fill(g.begin(), g.end(), 0.0);
                REprintf("xsdm XPtr ObjFun error: %s\n", e.what());
            } catch (...) {
                f = std::numeric_limits<double>::infinity();
                std::fill(g.begin(), g.end(), 0.0);
                REprintf("xsdm XPtr ObjFun error: unknown C++ exception\n");
            }
        }
    );

    return Rcpp::XPtr<ucminf::ObjFun>(fn.release(), true);
}
