#' @keywords internal
#' @noRd
NULL

# ---------------------------------------------------------------------------
# Pure-R reference implementations of the log-likelihood chain.
#
# These are kept for parity testing only -- the canonical (exported) names
# `loglik_math`, `loglik_bio`, `log_prob_detect` are thin R wrappers around
# the C++ implementations in `src/`. The `_r` references below were the
# pre-port pure-R code; they are not exported and should never be called
# from the optimizer hot path. They exist so that
# `tests/testthat/test-*_r_vs_cpp.R` can assert numerical equality
# between the R reference and the C++ wrapper.
# ---------------------------------------------------------------------------

#' Pure-R reference for `log_prob_detect`
#'
#' @keywords internal
#' @noRd
log_prob_detect_r <- function(env_dat,
                              mu,
                              sigltil,
                              sigrtil,
                              o_mat,
                              ctil,
                              pd,
                              return_prob = FALSE,
                              num_threads = RcppParallel::defaultNumThreads()) {
  check_env_array(env_dat)
  old <- RcppParallel::defaultNumThreads()
  on.exit(RcppParallel::setThreadOptions(numThreads = old), add = TRUE)
  RcppParallel::setThreadOptions(numThreads = num_threads)

  h <- like_neg_ltsgr(
    env_dat = env_dat,
    mu = mu,
    sigltil = sigltil,
    sigrtil = sigrtil,
    o_mat = o_mat
  )

  logpdetect <- log(pd) - log1pexp(ctil + h)

  if (return_prob) {
    logpdetect <- exp(logpdetect)
  }
  logpdetect
}

#' Pure-R reference for `loglik_bio`
#'
#' @keywords internal
#' @noRd
loglik_bio_r <- function(env_dat,
                         occ,
                         mu,
                         sigltil,
                         sigrtil,
                         o_mat,
                         ctil,
                         pd,
                         return_prob = FALSE,
                         sum_log_p = TRUE,
                         num_threads = RcppParallel::defaultNumThreads()) {
  checkmate::assert(
    checkmate::check_logical(occ, any.missing = FALSE),
    checkmate::check_integerish(occ, lower = 0, upper = 1, any.missing = FALSE),
    .var.name = "occ"
  )
  check_env_array(env_dat)
  checkmate::assert_numeric(mu, any.missing = FALSE, min.len = 1)
  checkmate::assert_numeric(sigltil, any.missing = FALSE, min.len = 1)
  checkmate::assert_numeric(sigrtil, any.missing = FALSE, min.len = 1)
  checkmate::assert_numeric(ctil, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(pd, any.missing = FALSE, len = 1)
  checkmate::assert_matrix(
    o_mat,
    min.rows = 1, min.cols = 1, any.missing = FALSE
  )

  RcppParallel::setThreadOptions(numThreads = num_threads)
  on.exit(
    RcppParallel::setThreadOptions(numThreads = RcppParallel::defaultNumThreads()),
    add = TRUE
  )

  log_p <- log_prob_detect_r(
    env_dat = env_dat,
    mu = mu,
    sigltil = sigltil,
    sigrtil = sigrtil,
    o_mat = o_mat,
    ctil = ctil,
    pd = pd,
    return_prob = FALSE
  )

  if (sum_log_p) {
    res <- sum(occ * log_p + (1 - occ) * log1mexp(-log_p))
  } else {
    res <- occ * log_p + (1 - occ) * log1mexp(-log_p)
  }

  if (return_prob) {
    res <- exp(res)
  }

  res
}

#' Pure-R reference for `loglik_math`
#'
#' @keywords internal
#' @noRd
loglik_math_r <- function(param_vector,
                          env_dat,
                          occ,
                          mask = NULL,
                          num_threads = RcppParallel::defaultNumThreads(),
                          negative = TRUE) {
  checkmate::assert(
    checkmate::check_logical(occ, any.missing = FALSE),
    checkmate::check_integerish(occ, lower = 0, upper = 1, any.missing = FALSE),
    .var.name = "occ"
  )
  check_env_array(env_dat)
  param_vector <- unlist(param_vector)
  checkmate::assert_vector(param_vector, any.missing = FALSE)
  p <- dim(env_dat)[3]

  if (is.null(names(param_vector))) {
    all_canonical <- names(make_mask_names(p))
    free_names <- if (!is.null(mask)) setdiff(all_canonical, names(mask)) else all_canonical
    if (length(param_vector) == length(free_names)) {
      names(param_vector) <- free_names
    } else {
      stop(
        "`param_vector` is unnamed and its length (", length(param_vector), ") ",
        "does not match the number of free parameters (", length(free_names), "). ",
        "Expected free parameters: ", paste(free_names, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  param_vector <- create_param_vector_masked_r(param_vector, mask, p)

  param_list <- math_to_bio_r(param_vector)

  checkmate::assert_numeric(param_list$mu, any.missing = FALSE, min.len = 1)
  checkmate::assert_numeric(param_list$sigltil, any.missing = FALSE, min.len = 1)
  checkmate::assert_numeric(param_list$sigrtil, any.missing = FALSE, min.len = 1)
  checkmate::assert_numeric(param_list$ctil, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(param_list$pd, any.missing = FALSE, len = 1)

  res <- loglik_bio_r(
    env_dat = env_dat,
    occ = occ,
    mu = param_list$mu,
    sigltil = param_list$sigltil,
    sigrtil = param_list$sigrtil,
    ctil = param_list$ctil,
    pd = param_list$pd,
    o_mat = param_list$o_mat,
    num_threads = num_threads
  )

  if (!negative) res else -res
}

# ---------------------------------------------------------------------------
# Pure-R reference implementations for the math <-> bio scale machinery.
#
# `math_to_bio_r` and `create_param_vector_masked_r` are the pre-port pure-R
# implementations. They are non-exported and exist solely as references for
# the parity tests in tests/testthat/test-*_r_vs_cpp.R, which assert that
# the C++ implementations (.math_to_bio_cpp, .build_canonical_param_vector_cpp)
# produce identical results.
# ---------------------------------------------------------------------------

#' Pure-R reference for `math_to_bio`
#'
#' @keywords internal
#' @noRd
math_to_bio_r <- function(param_vector) {
  # ---- Validate input structure ----
  checkmate::assert_numeric(param_vector,
                            names = "named",
                            any.missing = FALSE,
                            min.len = 5  # smallest possible length (p = 1)
  )

  # Infer p from the length
  n <- length(param_vector)
  p <- num_env_var(n)  # will error if n is not a valid num_par(p)

  # Verify that names exactly match the canonical order
  expected_names <- names(make_mask_names(p))
  if (!identical(names(param_vector), expected_names)) {
    stop(
      "`param_vector` names do not match the canonical order for p = ", p, ".\n",
      "Expected: ", paste(expected_names, collapse = ", "), "\n",
      "Received: ", paste(names(param_vector), collapse = ", ")
    )
  }
  param_vector <- unlist(param_vector)
  # ---- Extract components using canonical names ----
  mu       <- param_vector[grep("^mu[0-9]+$", names(param_vector))] |> as.numeric()
  sigltil  <- param_vector[grep("^sigltil", names(param_vector))]   |> exp() |> as.numeric()
  sigrtil  <- param_vector[grep("^sigrtil", names(param_vector))]   |> exp() |> as.numeric()
  ctil     <- param_vector[grep("^ctil", names(param_vector))]      |> as.numeric()
  pd       <- param_vector[grep("^pd", names(param_vector))]        |> as.numeric() |> expit()
  o_par    <- param_vector[grep("^o_par", names(param_vector))]     |> as.numeric()

  o_mat <- build_orthogonal_matrix(if (length(o_par) == 0) NULL else o_par)

  list(
    mu      = mu,
    sigltil = sigltil,
    sigrtil = sigrtil,
    ctil    = ctil,
    pd      = pd,
    o_mat   = o_mat
  )
}

#' Pure-R reference for `create_param_vector_masked`
#'
#' @keywords internal
#' @noRd
create_param_vector_masked_r <- function(param_vector, mask = NULL, p) {
  checkmate::assert_count(p, positive = TRUE)

  # Build full canonical skeleton (with mask applied first)
  out <- create_mask(mask = mask, p = p)
  allowed_names <- names(out)
  if (is.null(allowed_names)) {
    stop("`create_mask(mask = NULL, p)` must return a named vector;
         names are canonical.")
  }

  # Validate param_vector (required)
  if (is.null(param_vector)) {
    stop("`param_vector` must not be NULL.")
  }
  checkmate::assert_numeric(param_vector, any.missing = FALSE, names = "named")
  bad_pv <- setdiff(names(param_vector), allowed_names)
  if (length(bad_pv) > 0) {
    stop(
      "Unexpected name(s) in `param_vector`: ", paste(bad_pv, collapse = ", "),
      ". Allowed: ", paste(allowed_names, collapse = ", "), "."
    )
  }

  # Validate mask when present and enforce disjointness
  if (!is.null(mask)) {
    checkmate::assert_numeric(mask, any.missing = FALSE, names = "named")
    bad_mask <- setdiff(names(mask), allowed_names)
    if (length(bad_mask) > 0) {
      stop(
        "Unexpected name(s) in `mask`: ", paste(bad_mask, collapse = ", "),
        ". Allowed: ", paste(allowed_names, collapse = ", "), "."
      )
    }
    overlap <- intersect(names(param_vector), names(mask))
    if (length(overlap) > 0) {
      stop(
        "`param_vector` and `mask` must be complementary (disjoint). ",
        "Overlapping names: ", paste(overlap, collapse = ", "), "."
      )
    }
  }

  out[names(param_vector)] <- param_vector

  if (anyNA(out)) {
    missing_names <- names(out)[is.na(out)]
    stop(
      "The final parameter vector contains NA values for: ",
      paste(missing_names, collapse = ", "),
      ". Provide values via `mask` and/or `param_vector` so all positions are
      filled."
    )
  }

  out
}

# ---------------------------------------------------------------------------
# Small non-exported helpers consolidated from individual R/*.R files.
#
# These were previously each in their own file (R/check_env_array.R,
# R/logit.R, R/permutations.R). They are kept internal-only and grouped
# here for ease of maintenance.
# ---------------------------------------------------------------------------

#' Validate the environmental data array
#'
#' Checks that \code{env_dat} is a 3-dimensional array with dimensions
#' \code{n_loc x n_time x p} (number of locations x time-series length x
#' number of environmental variables) and contains no missing values.
#' Throws an informative error if any condition is violated.
#'
#' @keywords internal
#' @noRd
check_env_array <- function(env_dat, name = "env_dat") {
  checkmate::assert_array(
    env_dat,
    d = 3,
    any.missing = FALSE,
    .var.name = name
  )
  invisible(env_dat)
}

#' Logit function. Inverse of the expit function.
#'
#' @keywords internal
#' @noRd
logit <- function(x) {
  no <- (x < 0) | (x > 1)
  out <- numeric(length(x))
  out[no] <- NaN
  out[!no] <- log(x[!no] / (1 - x[!no]))
  dim(out) <- dim(x)
  out
}

#' Enumerate Permutations of Vector Elements
#'
#' Adapted from \pkg{gtools::permutations} (Bill Venables; extended by
#' Gregory R. Warnes to handle \code{repeats.allowed}). Used internally
#' by \code{distance_between_params} for the orthogonal-matrix
#' equivalence-class search.
#'
#' @keywords internal
#' @noRd
permutations <- function(n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE) {
  if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n %% 1) !=
    0) {
    stop("bad value of n")
  }
  if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r %% 1) !=
    0) {
    stop("bad value of r")
  }
  if (!is.atomic(v) || length(v) < n) {
    stop("v is either non-atomic or too short")
  }
  if ((r > n) & repeats.allowed == FALSE) {
    stop("r > n and repeats.allowed=FALSE")
  }
  if (set) {
    v <- unique(sort(v))
    if (length(v) < n) {
      stop("too few different elements")
    }
  }
  v0 <- vector(mode(v), 0)
  if (repeats.allowed) {
    sub <- function(n, r, v) {
      if (r == 1) {
        matrix(v, n, 1)
      } else if (n == 1) {
        matrix(v, 1, r)
      } else {
        inner <- Recall(n, r - 1, v)
        cbind(rep(v, rep(nrow(inner), n)), matrix(t(inner),
          ncol = ncol(inner), nrow = nrow(inner) * n,
          byrow = TRUE
        ))
      }
    }
  } else {
    sub <- function(n, r, v) {
      if (r == 1) {
        matrix(v, n, 1)
      } else if (n == 1) {
        matrix(v, 1, r)
      } else {
        X <- NULL
        for (i in 1:n) {
          X <- rbind(X, cbind(v[i], Recall(n -
            1, r - 1, v[-i])))
        }
        X
      }
    }
  }
  sub(n, r, v[1:n])
}
