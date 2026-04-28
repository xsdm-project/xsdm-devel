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

  h <- like_neg_ltsgr_cpp(
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

  param_vector <- create_param_vector_masked(param_vector, mask, p)

  param_list <- math_to_bio(param_vector)

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
