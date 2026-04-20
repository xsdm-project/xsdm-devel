//' Compute likelihood for LTSG model
//'
//' This function calculates a likelihood-like measure using orthogonal matrices,
//' environmental data, and diagonal matrices, leveraging parallel computation.
//'
//' @name like_ltsg
//' @title Compute likelihood for LTSG model
//' @param mu Numeric vector of means (length equal to number of rows in `env_m`)
//' @param env_m Numeric matrix of environmental data
//' @param dl_mat Diagonal matrix (as NumericMatrix)
//' @param drl_mat Diagonal matrix (as NumericMatrix)
//' @param ortho_m Numeric matrix (orthogonal basis)
//' @param q Integer, number of rows for reshaping
//' @param r Integer, number of columns for reshaping.
//' @return A numeric vector of length `r` with computed sums.

#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

// Worker struct for parallel computation over columns of env_m
struct ColumnWorker : public Worker {
   // Views into R objects for efficient access
   const RMatrix<double> ortho_m;   // Orthogonal matrix (k x p)
   const RMatrix<double> env_m;     // Environmental data (p x T)
   const RMatrix<double> dl_mat;    // Diagonal scaling matrix (k x k)
   const RMatrix<double> drl_mat;   // Diagonal scaling matrix for positive part
   const RVector<double> mu;        // Mean vector (length p)
   RVector<double> output;          // Output vector (length T)
   
   // Constructor: bind R objects to worker
   ColumnWorker(const NumericMatrix& ortho_m,
                const NumericMatrix& env_m,
                const NumericMatrix& dl_mat,
                const NumericMatrix& drl_mat,
                NumericVector& mu,
                NumericVector& output)
     : ortho_m(ortho_m), env_m(env_m), dl_mat(dl_mat), drl_mat(drl_mat), mu(mu),
       output(output) {}
   
   // Parallel operator: compute for columns in [begin, end)
   void operator()(std::size_t begin, std::size_t end) {
     int nrow_ortho_m = ortho_m.nrow(); // Number of orthogonal directions (k)
     int ncol_ortho_m = ortho_m.ncol(); // Number of variables (p)
     
     for (std::size_t j = begin; j < end; j++) {
       double col_sum = 0.0; // Accumulator for column j
       
       // Loop over each orthogonal row i
       for (int i = 0; i < nrow_ortho_m; i++) {
         double dot_product = 0.0;
         
         // Compute projection: sum over variables
         for (int k = 0; k < ncol_ortho_m; k++) {
           dot_product += ortho_m(i, k) * (env_m(k, j) - mu[k]);
         }
         
         // Apply asymmetric scaling using dl_mat and drl_mat
         double usym = (dl_mat(i, i) * dot_product +
                        drl_mat(i, i) * std::max(0.0, dot_product));
         
         // Accumulate squared contribution
         col_sum += usym * usym;
       }
       
       // Store result for column j
       output[j] = col_sum;
     }
   }
 };
 
// [[Rcpp::export]]
NumericVector like_ltsg(NumericVector mu,
                     NumericMatrix env_m,
                     NumericMatrix dl_mat,
                     NumericMatrix drl_mat,
                     NumericMatrix ortho_m,
                     int q,
                     int r) {
 // Check matrix compatibility: ortho_m columns must match env_m rows
 if (ortho_m.ncol() != env_m.nrow()) {
   stop("Matrix dimensions are not compatible for multiplication.");
 }
 
 int ncol_env_m = env_m.ncol(); // Number of observations (T)
 
 // Assumes env_m is column-major with time varying fastest:
 // column j corresponds to (location, time) via idx = j*q + i
 
 if (q * r != ncol_env_m)
   stop("q * r must equal env_m.ncol().");
 NumericVector output(ncol_env_m); // Initialize output vector
 
 // Create worker and run parallel computation
 ColumnWorker worker(ortho_m, env_m, dl_mat, drl_mat, mu, output);
 parallelFor(0, ncol_env_m, worker);
 
 // Step 2: Reshape output into q x r matrix
 NumericMatrix d_mat(q, r);
 for (int j = 0; j < r; j++) {
   for (int i = 0; i < q; i++) {
     int idx = j * q + i; // Compute index in output
     d_mat(i, j) = output[idx] / (2 * q); // Normalize by 2*q
   }
 }
   
 // Step 3: Sum columns of reshaped matrix to get final result
 NumericVector result(r);
 for (int j = 0; j < r; j++) {
   double sum = 0.0;
   for (int i = 0; i < q; i++) {
     sum += d_mat(i, j);
   }
   result[j] = sum;
 }
   
 return result; // Return vector of length r
}

