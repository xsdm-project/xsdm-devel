# inst/benchmarks/loglik_cpp_vs_r.R
#
# Timing comparison: pure-C++ loglik_math_cpp vs the R reference
# loglik_math. The XPtr hot path calls loglik_math_eval (which
# loglik_math_cpp wraps), so this is a proxy for the per-iteration
# speedup obtained by make_loglik_math_xptr_cpp over the legacy
# R-callback ObjFun.
#
# Usage:
#   Rscript inst/benchmarks/loglik_cpp_vs_r.R
#
# Uses base R's system.time so it runs without extra dependencies.

suppressPackageStartupMessages({
  library(xsdm)
})

one_case <- function(n, Tt, p, reps = 50L) {
  set.seed(7L)
  env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
  occ     <- rep(c(1L, 0L), length.out = n)

  pv <- make_mask_names(p)
  pv[] <- rnorm(length(pv), sd = 0.3)

  t_r <- system.time({
    for (i in seq_len(reps)) {
      loglik_math(pv, env_dat, occ, negative = TRUE,
                  num_threads = 1L)
    }
  })[["elapsed"]]

  t_cpp <- system.time({
    for (i in seq_len(reps)) {
      loglik_math_cpp(pv, env_dat, occ, negative = TRUE,
                      num_threads = 1L)
    }
  })[["elapsed"]]

  data.frame(
    case    = sprintf("n=%d,Tt=%d,p=%d", n, Tt, p),
    reps    = reps,
    r_sec   = round(t_r,   4),
    cpp_sec = round(t_cpp, 4),
    speedup = round(t_r / t_cpp, 2)
  )
}

rows <- rbind(
  one_case(  10,   5, 2, reps = 200L),
  one_case(  50,  10, 2, reps = 100L),
  one_case( 200,  20, 2, reps =  50L),
  one_case(1000,  20, 2, reps =  10L)
)

cat("loglik_math_cpp (pure C++) vs loglik_math (R reference)\n")
cat(strrep("=", 60), "\n", sep = "")
print(rows, row.names = FALSE)
