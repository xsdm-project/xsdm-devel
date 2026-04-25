// src/expm_skew.hpp
//
// Matrix exponential of a real skew-symmetric p x p matrix S = L - L^T,
// where L is the strictly-lower-triangular matrix whose entries are taken
// from a vector in column-major order.
//
// Mirrors R/build_orthogonal_matrix.R:
//
//   s_matrix[lower.tri(s_matrix)] <- entries   # column-major in R
//   s_matrix <- s_matrix - t(s_matrix)
//   o_matrix <- expm::expm(s_matrix)
//
// Algorithm: scaling-and-squaring Padé-13, following Higham (2005),
// "The Scaling and Squaring Method for the Matrix Exponential Revisited".
// For the small p that xsdm ships with (p <= ~10 in practice) the cost is
// negligible next to the per-cell likelihood evaluation, and the relative
// accuracy is ~1e-14.
//
// For p == 2 there is a closed-form rotation, used as a fast path.
// For p == 1 the matrix is 1x1 identity.
//
// This header is header-only so both the R-callable wrapper
// (math_to_bio_cpp) and the XPtr ObjFun closure can include it without
// a separate .cpp file.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace xsdm {

// ---------------------------------------------------------------------------
// Small dense column-major BLAS-lite helpers.
//
// A p x p matrix A is stored as std::vector<double>(p*p) with
//   A(i, j) = A_data[i + j * p]      (column-major, matches R / LAPACK)
// ---------------------------------------------------------------------------

// C = A * B   (p x p times p x p)
inline void mm_cm(const double* A, const double* B, double* C, int p) {
    for (int j = 0; j < p; ++j) {
        for (int i = 0; i < p; ++i) {
            double s = 0.0;
            for (int k = 0; k < p; ++k) {
                s += A[i + k * p] * B[k + j * p];
            }
            C[i + j * p] = s;
        }
    }
}

// A <- A + alpha * B
inline void axpy_mat(double* A, const double* B, double alpha, int p) {
    const int n = p * p;
    for (int i = 0; i < n; ++i) A[i] += alpha * B[i];
}

// A <- scalar * I   (overwrite)
inline void set_scaled_identity(double* A, double scalar, int p) {
    const int n = p * p;
    std::fill(A, A + n, 0.0);
    for (int i = 0; i < p; ++i) A[i + i * p] = scalar;
}

// inf-norm = max row absolute-row-sum.
inline double inf_norm(const double* A, int p) {
    double best = 0.0;
    for (int i = 0; i < p; ++i) {
        double row = 0.0;
        for (int j = 0; j < p; ++j) row += std::abs(A[i + j * p]);
        best = std::max(best, row);
    }
    return best;
}

// Solve (I - U/2 + V/2) X = (I + U/2 + V/2) type system is simpler if we
// form P = V + U, Q = V - U  (see Higham 2005 eq. (2.1)); actually the
// classic formulation for [m/m]-Padé of exp is:
//     exp(A) ≈ Q(A)^{-1} P(A)
// where  P(A) = U(A) + V(A),  Q(A) = -U(A) + V(A)
// and U, V are defined below (Padé-13).
//
// We implement Gauss-Jordan elimination with partial pivoting on the p x p
// system Q * X = P to keep the file dependency-free (no LAPACK).

// In-place solve Q * X = P, overwriting X with the solution.
// Returns false if Q is numerically singular.
inline bool gauss_solve(double* Q, double* X, int p) {
    // augment [Q | X] and eliminate column by column.
    std::vector<double> M(static_cast<std::size_t>(p) * 2 * p, 0.0);
    for (int j = 0; j < p; ++j) {
        for (int i = 0; i < p; ++i) {
            M[i + j * p]            = Q[i + j * p];
            M[i + (j + p) * p]      = X[i + j * p];
        }
    }

    for (int k = 0; k < p; ++k) {
        // Partial pivot: find row with max |M[r, k]| for r >= k.
        int piv = k;
        double best = std::abs(M[k + k * p]);
        for (int r = k + 1; r < p; ++r) {
            double v = std::abs(M[r + k * p]);
            if (v > best) { best = v; piv = r; }
        }
        if (best < 1e-300) return false;
        if (piv != k) {
            for (int j = 0; j < 2 * p; ++j) {
                std::swap(M[k + j * p], M[piv + j * p]);
            }
        }
        // Normalize pivot row.
        double inv = 1.0 / M[k + k * p];
        for (int j = k; j < 2 * p; ++j) M[k + j * p] *= inv;
        // Eliminate other rows.
        for (int r = 0; r < p; ++r) {
            if (r == k) continue;
            double f = M[r + k * p];
            if (f == 0.0) continue;
            for (int j = k; j < 2 * p; ++j) M[r + j * p] -= f * M[k + j * p];
        }
    }

    for (int j = 0; j < p; ++j) {
        for (int i = 0; i < p; ++i) {
            X[i + j * p] = M[i + (j + p) * p];
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// expm_pade13: matrix exponential of a small dense column-major p x p matrix.
// ---------------------------------------------------------------------------
inline void expm_pade13(const double* A_in, double* E_out, int p) {
    const int n = p * p;

    // Padé-13 coefficients (Higham 2005, Table 2.3).
    static const double b[] = {
        64764752532480000.0,
        32382376266240000.0,
         7771770303897600.0,
         1187353796428800.0,
          129060195264000.0,
           10559470521600.0,
             670442572800.0,
              33522128640.0,
               1323241920.0,
                 40840800.0,
                   960960.0,
                    16380.0,
                      182.0,
                        1.0
    };

    // Scaling: choose s so that ||A / 2^s||_inf <= theta_13.
    const double theta13 = 5.371920351148152;  // Higham 2005
    std::vector<double> A(A_in, A_in + n);
    double norm = inf_norm(A.data(), p);
    int s = 0;
    if (norm > theta13) {
        s = static_cast<int>(std::ceil(std::log2(norm / theta13)));
        if (s < 0) s = 0;
        double scale = std::ldexp(1.0, -s);   // 2^-s
        for (int i = 0; i < n; ++i) A[i] *= scale;
    }

    // Compute A^2, A^4, A^6.
    std::vector<double> A2(n), A4(n), A6(n);
    mm_cm(A.data(),  A.data(),  A2.data(), p);
    mm_cm(A2.data(), A2.data(), A4.data(), p);
    mm_cm(A2.data(), A4.data(), A6.data(), p);

    // U = A * (A^6 * (b13 A^6 + b11 A^4 + b9 A^2) + b7 A^6 + b5 A^4 + b3 A^2 + b1 I)
    // V =     A^6 * (b12 A^6 + b10 A^4 + b8 A^2) + b6 A^6 + b4 A^4 + b2 A^2 + b0 I
    std::vector<double> tmp(n), U(n, 0.0), V(n, 0.0);

    // inner-U = b13 A^6 + b11 A^4 + b9 A^2
    std::vector<double> innerU(n, 0.0);
    axpy_mat(innerU.data(), A6.data(), b[13], p);
    axpy_mat(innerU.data(), A4.data(), b[11], p);
    axpy_mat(innerU.data(), A2.data(), b[9],  p);

    // tmp = A^6 * inner-U
    mm_cm(A6.data(), innerU.data(), tmp.data(), p);

    // U-rhs = tmp + b7 A^6 + b5 A^4 + b3 A^2 + b1 I
    std::vector<double> Urhs(tmp);
    axpy_mat(Urhs.data(), A6.data(), b[7], p);
    axpy_mat(Urhs.data(), A4.data(), b[5], p);
    axpy_mat(Urhs.data(), A2.data(), b[3], p);
    for (int i = 0; i < p; ++i) Urhs[i + i * p] += b[1];

    // U = A * Urhs
    mm_cm(A.data(), Urhs.data(), U.data(), p);

    // inner-V = b12 A^6 + b10 A^4 + b8 A^2
    std::vector<double> innerV(n, 0.0);
    axpy_mat(innerV.data(), A6.data(), b[12], p);
    axpy_mat(innerV.data(), A4.data(), b[10], p);
    axpy_mat(innerV.data(), A2.data(), b[8],  p);

    // V = A^6 * inner-V + b6 A^6 + b4 A^4 + b2 A^2 + b0 I
    mm_cm(A6.data(), innerV.data(), V.data(), p);
    axpy_mat(V.data(), A6.data(), b[6], p);
    axpy_mat(V.data(), A4.data(), b[4], p);
    axpy_mat(V.data(), A2.data(), b[2], p);
    for (int i = 0; i < p; ++i) V[i + i * p] += b[0];

    // P = V + U,  Q = V - U,  solve Q * R = P for R = Padé-13(A).
    std::vector<double> P(V), Q(V);
    axpy_mat(P.data(), U.data(),  1.0, p);
    axpy_mat(Q.data(), U.data(), -1.0, p);

    std::vector<double> R(P);  // will be overwritten by gauss_solve
    gauss_solve(Q.data(), R.data(), p);

    // Undo the scaling: R <- R^{2^s}
    for (int k = 0; k < s; ++k) {
        mm_cm(R.data(), R.data(), tmp.data(), p);
        R.swap(tmp);
    }

    std::copy(R.begin(), R.end(), E_out);
}

// ---------------------------------------------------------------------------
// build_orthogonal_matrix_cpp
//
// Identical semantics to R/build_orthogonal_matrix.R:
//   - entries has length q = p * (p - 1) / 2, in column-major lower-triangular
//     order (i.e., the same order that R's `s_matrix[lower.tri(s_matrix)] <-
//     entries` produces).
//   - Build S (skew-symmetric) = L - L^T, then O = expm(S).
//
// For p == 1:  returns 1, entries may be nullptr (q == 0).
// For p == 2:  closed-form rotation [[cos θ, sin θ], [-sin θ, cos θ]] where
//              θ = entries[0].  Relative error ~ double machine eps.
// For p >= 3:  Padé-13 scaling-and-squaring.
// ---------------------------------------------------------------------------
inline void build_orthogonal_matrix_cpp(
    const double* entries,
    int           p,
    double*       o_mat_out      // length p*p, column-major
) {
    if (p == 1) {
        o_mat_out[0] = 1.0;
        return;
    }
    if (p == 2) {
        const double theta = entries[0];
        const double c = std::cos(theta);
        const double s = std::sin(theta);
        // Match R: L has entries[0] at (2,1); S = L - L^T so
        //   S = [[0, -theta], [theta, 0]]
        // => expm(S) = [[cos θ, -sin θ], [sin θ, cos θ]] (column-major below)
        o_mat_out[0] =  c;   // (1,1)
        o_mat_out[1] =  s;   // (2,1)
        o_mat_out[2] = -s;   // (1,2)
        o_mat_out[3] =  c;   // (2,2)
        return;
    }

    // General case: build S in column-major, then Padé-13.
    const int n = p * p;
    std::vector<double> S(n, 0.0);

    // Fill strict lower triangle in column-major order, matching R:
    //   lower.tri() visits (row=2..p, col=1), (row=3..p, col=2), ...
    int k = 0;
    for (int j = 0; j < p; ++j) {
        for (int i = j + 1; i < p; ++i) {
            S[i + j * p] =  entries[k];
            S[j + i * p] = -entries[k];
            ++k;
        }
    }

    expm_pade13(S.data(), o_mat_out, p);
}

} // namespace xsdm
