// src/dist_between_params.cpp
//
// Pure-C++ backend for R/dist_between_params.R. The R wrapper continues to
// perform all input validation and math->biological-scale conversion; this
// file is responsible for the two numerically heavy pieces:
//
//   (1) Build the (dd x dd) cost matrix and the +/-1 sign matrix that
//       records which orientation of each pairing (pos vs neg, see R doc)
//       was chosen.
//   (2) Solve the resulting linear sum assignment problem (LSAP) via the
//       classical O(n^3) Hungarian algorithm, returning the optimal
//       column-to-row permutation.
//
// -----------------------------------------------------------------------------
// Algorithmic credit and licence note
// -----------------------------------------------------------------------------
// The linear sum assignment algorithm implemented in `solve_lsap_core` below
// is the classical Hungarian / Kuhn-Munkres algorithm, described in:
//
//   * H. W. Kuhn (1955). "The Hungarian Method for the Assignment Problem."
//     Naval Research Logistics Quarterly, 2(1-2), 83-97.
//   * J. Munkres (1957). "Algorithms for the Assignment and Transportation
//     Problems." Journal of the SIAM, 5(1), 32-38.
//   * R. Jonker, A. Volgenant (1987). "A Shortest Augmenting Path Algorithm
//     for Dense and Sparse Linear Assignment Problems." Computing 38, 325-340.
//
// The implementation below is a clean-room write-up of the potentials-based
// O(n^3) variant described in those papers. It is *not* derived from the
// source code of any existing package. In particular, we do not reuse code
// from `clue::solve_LSAP` (Hornik, Boehm; CRAN package `clue`, GPL-2), even
// though we match its numerical output in the parity tests. Keeping this
// implementation independent avoids the GPL-2 / AGPL-3 licence mismatch and
// lets the package remain under its declared AGPL (>= 3).
//
// Users interested in the original R implementation of the same LSAP
// algorithm are encouraged to look at the `clue` package on CRAN.
// -----------------------------------------------------------------------------

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

// Solve the (square) linear sum assignment problem on a n x n cost matrix
// laid out row-major: cost[i * n + j]. Returns a vector `assign` of length
// n such that row i is matched to column assign[i] (0-indexed internally).
//
// Implementation: potentials-based O(n^3) Hungarian. Uses +1 dummy source
// row (index 0); rows/cols below are 1-indexed inside this routine.
std::vector<int> solve_lsap_core(const std::vector<double>& cost, int n) {
    const double INF = std::numeric_limits<double>::infinity();

    // Potentials u (rows), v (columns). Size (n+1) so we can 1-index.
    std::vector<double> u(n + 1, 0.0);
    std::vector<double> v(n + 1, 0.0);

    // p[j] = row (1..n) currently assigned to column j; 0 means unassigned.
    std::vector<int> p(n + 1, 0);

    // way[j] = predecessor column in the alternating path, for reconstruction.
    std::vector<int> way(n + 1, 0);

    for (int i = 1; i <= n; ++i) {
        // Start augmenting-path search rooted at row i (stored in p[0]).
        p[0] = i;
        int j0 = 0;

        std::vector<double> minv(n + 1, INF);
        std::vector<char>   used(n + 1, false);

        do {
            used[j0] = true;
            const int i0 = p[j0];
            double delta = INF;
            int j1 = -1;

            for (int j = 1; j <= n; ++j) {
                if (used[j]) continue;
                // Reduced cost using current potentials.
                const double cur = cost[(i0 - 1) * n + (j - 1)] - u[i0] - v[j];
                if (cur < minv[j]) {
                    minv[j] = cur;
                    way[j]  = j0;
                }
                if (minv[j] < delta) {
                    delta = minv[j];
                    j1 = j;
                }
            }

            // Update potentials along the explored fragment.
            for (int j = 0; j <= n; ++j) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j]    -= delta;
                } else {
                    minv[j] -= delta;
                }
            }

            j0 = j1;
        } while (p[j0] != 0);

        // Reverse the alternating path to commit the matching.
        while (j0 != 0) {
            const int j1 = way[j0];
            p[j0] = p[j1];
            j0    = j1;
        }
    }

    // Translate column->row assignment into row->column assignment.
    std::vector<int> assign(n, -1);
    for (int j = 1; j <= n; ++j) {
        if (p[j] >= 1 && p[j] <= n) assign[p[j] - 1] = j - 1;
    }
    return assign;
}

// Inverse-scale squared difference used for sigltil / sigrtil comparisons.
// Matches the R helper `sigdistsq <- function(x, y) (1/x - 1/y)^2` including
// IEEE semantics for sigltil/sigrtil entries of +Inf (1/Inf == 0).
inline double sigdistsq(double x, double y) {
    const double ix = std::isinf(x) ? 0.0 : 1.0 / x;
    const double iy = std::isinf(y) ? 0.0 : 1.0 / y;
    const double d  = ix - iy;
    return d * d;
}

} // anonymous namespace

// ----------------------------------------------------------------------------
// .solve_lsap_cpp
//
// Standalone Rcpp wrapper over the Hungarian solver. Exposed with a dot
// prefix (unexported). Primarily useful for parity testing against
// clue::solve_LSAP; the main dist_between_params wrapper calls the
// C++-side implementation directly and does not go through this.
// ----------------------------------------------------------------------------

//' @keywords internal
//' @noRd
// [[Rcpp::export(.solve_lsap_cpp)]]
Rcpp::IntegerVector solve_lsap_cpp(Rcpp::NumericMatrix cost) {
    const int n = cost.nrow();
    if (cost.ncol() != n) {
        Rcpp::stop("solve_lsap_cpp: `cost` must be a square matrix.");
    }
    if (n == 0) {
        return Rcpp::IntegerVector::create();
    }

    // Flatten to row-major for the solver.
    std::vector<double> flat(static_cast<std::size_t>(n) * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const double v = cost(i, j);
            if (!R_FINITE(v)) {
                Rcpp::stop("solve_lsap_cpp: `cost` must contain only finite values.");
            }
            flat[static_cast<std::size_t>(i) * n + j] = v;
        }
    }

    const std::vector<int> assign = solve_lsap_core(flat, n);
    Rcpp::IntegerVector out(n);
    for (int i = 0; i < n; ++i) {
        // Return 1-indexed to match clue::solve_LSAP output convention.
        out[i] = assign[i] + 1;
    }
    return out;
}

// ----------------------------------------------------------------------------
// .dist_between_params_cpp
//
// Build the cost + posneg matrices from the biological-scale sigltil / sigrtil
// / o_mat pairs of two parameter sets, solve the LSAP, and return a list with
//
//   * distance   : sqrt( sum(cost at pairing) + sum((mu1-mu2)^2)
//                        + (ctil1-ctil2)^2 + (pd1-pd2)^2 )
//   * perm       : integer vector, 1-indexed; perm[i] is the column of p1
//                  matched to column i of p2
//   * posneg     : integer matrix (dd x dd), entries +1 or -1 recording which
//                  orientation (pos / neg) was chosen for each candidate pair;
//                  ties are broken to +1 to match the R reference exactly
//   * cost       : double matrix (dd x dd), the assignment-problem cost matrix
//                  (kept for diagnostic / parity tests)
//
// All arguments are biological-scale. No validation beyond dimension checks
// is performed here; the R wrapper is responsible for user-facing checks.
// ----------------------------------------------------------------------------

//' @keywords internal
//' @noRd
// [[Rcpp::export(.dist_between_params_cpp)]]
Rcpp::List dist_between_params_cpp(
    Rcpp::NumericVector mu1,
    Rcpp::NumericVector sigltil1,
    Rcpp::NumericVector sigrtil1,
    Rcpp::NumericMatrix o_mat1,
    double ctil1,
    double pd1,
    Rcpp::NumericVector mu2,
    Rcpp::NumericVector sigltil2,
    Rcpp::NumericVector sigrtil2,
    Rcpp::NumericMatrix o_mat2,
    double ctil2,
    double pd2
) {
    const int dd = o_mat1.nrow();
    if (o_mat1.ncol() != dd || o_mat2.nrow() != dd || o_mat2.ncol() != dd) {
        Rcpp::stop("dist_between_params_cpp: o_mat1 and o_mat2 must be square "
                   "p x p matrices of the same dimension.");
    }
    if (mu1.size() != dd || mu2.size() != dd ||
        sigltil1.size() != dd || sigrtil1.size() != dd ||
        sigltil2.size() != dd || sigrtil2.size() != dd) {
        Rcpp::stop("dist_between_params_cpp: mu/sigltil/sigrtil vectors must "
                   "all have length dd (= nrow(o_mat1)).");
    }

    // Build cost + posneg matrices.
    Rcpp::NumericMatrix cost(dd, dd);
    Rcpp::IntegerMatrix posneg(dd, dd);
    std::vector<double> cost_flat(static_cast<std::size_t>(dd) * dd);

    for (int cc2 = 0; cc2 < dd; ++cc2) {
        for (int cc1 = 0; cc1 < dd; ++cc1) {
            // pos = sum_k (o2[k,cc2] - o1[k,cc1])^2
            //       + sigdistsq(sigltil2[cc2], sigltil1[cc1])
            //       + sigdistsq(sigrtil2[cc2], sigrtil1[cc1])
            double pos = 0.0;
            double neg = 0.0;
            for (int k = 0; k < dd; ++k) {
                const double a = o_mat2(k, cc2);
                const double b = o_mat1(k, cc1);
                const double dpos = a - b;
                const double dneg = a + b;
                pos += dpos * dpos;
                neg += dneg * dneg;
            }
            pos += sigdistsq(sigltil2[cc2], sigltil1[cc1]);
            pos += sigdistsq(sigrtil2[cc2], sigrtil1[cc1]);
            neg += sigdistsq(sigltil2[cc2], sigrtil1[cc1]);
            neg += sigdistsq(sigrtil2[cc2], sigltil1[cc1]);

            const double best = std::min(pos, neg);
            cost(cc2, cc1) = best;
            cost_flat[static_cast<std::size_t>(cc2) * dd + cc1] = best;

            // R reference: if pos < neg -> +1; if neg < pos -> -1; tie -> +1.
            posneg(cc2, cc1) = (neg < pos) ? -1 : 1;
        }
    }

    // Solve LSAP.
    const std::vector<int> assign_0 = solve_lsap_core(cost_flat, dd);
    Rcpp::IntegerVector perm(dd);
    for (int i = 0; i < dd; ++i) perm[i] = assign_0[i] + 1; // 1-indexed

    // Assemble distance.
    double sum_cost = 0.0;
    for (int i = 0; i < dd; ++i) {
        sum_cost += cost(i, assign_0[i]);
    }

    double sq_other = 0.0;
    for (int i = 0; i < dd; ++i) {
        const double d = mu1[i] - mu2[i];
        sq_other += d * d;
    }
    const double dc = ctil1 - ctil2;
    const double dp = pd1 - pd2;
    sq_other += dc * dc + dp * dp;

    const double distance = std::sqrt(sum_cost + sq_other);

    return Rcpp::List::create(
        Rcpp::Named("distance")             = distance,
        Rcpp::Named("perm")                 = perm,
        Rcpp::Named("posneg")               = posneg,
        Rcpp::Named("cost")                 = cost,
        Rcpp::Named("sq_dist_other_params") = sq_other
    );
}
