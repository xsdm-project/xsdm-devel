// src/solve_lsap.cpp
//
// Linear sum assignment problem (LSAP) solver used by R/dist_between_params.R.
// The R wrapper does all input validation, math->biological-scale conversion,
// and cost-matrix construction; this file is responsible for one piece only:
//
//   * Solve the resulting linear sum assignment problem on a square cost
//     matrix via the classical O(n^3) Hungarian algorithm, returning the
//     optimal column-to-row permutation.
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

} // anonymous namespace

// ----------------------------------------------------------------------------
// .solve_lsap_cpp
//
// Rcpp wrapper over the Hungarian solver. Exposed (unexported) under a dot
// prefix and called from R/dist_between_params.R as
// `as.integer(.solve_lsap_cpp(cost))`. Output is 1-indexed to match the
// `clue::solve_LSAP` convention so the two are interchangeable in R.
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
