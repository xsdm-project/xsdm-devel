#' Internal helper: run ucminf for one starting vector
#' @keywords internal
.optimize_loglik_math <- function(
    param_vector,  # starting values for FREE params (named)
    env_dat, occ, mask,   # data & fixed params
    num_threads,          # threads for loglik
    base_control,         # list passed to ucminf
    invh_lt = NULL,       # optional inverse Hessian LT warm-start
    optimizer_fun = ucminf::ucminf  # <- dependency injection for testing
) {
  # Merge invH warm-start if provided
  ctrl <- base_control
  if (!is.null(invh_lt)) ctrl$invhessian.lt <- invh_lt

  # Run optimizer with protection
  # NOTE: input coercion and validation are inside the tryCatch so that
  # non-finite starting values produce a structured ucminf_error instead of
  # an uncaught assertion failure (Fix 3).
  out <- tryCatch(
    {
      # Ensure numeric & finite par; keep names (inside tryCatch so that
      # non-finite starting values yield a structured ucminf_error — Fix 3)
      v <- as.numeric(param_vector)
      names(v) <- names(param_vector)
      checkmate::assert_numeric(v, any.missing = FALSE, finite = TRUE)
      param_vector <- v

      res <- optimizer_fun(
        par      = param_vector,
        fn       = loglik_math,
        env_dat  = env_dat,
        mask     = mask,
        occ      = occ,
        negative = TRUE,
        num_threads = num_threads,
        control  = ctrl,
        hessian  = FALSE
      )
      
      # Restore names if optimizer dropped them
      if (is.null(names(res$par)) && !is.null(names(param_vector))) {
        names(res$par) <- names(param_vector)
      }
      
      # If optimizer returned non-finite params, sanitize and mark custom code
      if (any(!is.finite(res$par))) {
        res$par[!is.finite(res$par)] <- param_vector[!is.finite(res$par)]
        res$convergence <- if (!is.null(res$convergence)) as.integer(res$convergence) else NA_integer_
        res$convergence <- -99L  # custom code for non-finite params
      }
      
      res
    },
    error = function(e) {
      # Structured failure: caller can still assemble results
      structure(
        list(
          par         = param_vector,
          value       = Inf,  # => loglik = -Inf
          convergence = NA_integer_,
          error       = conditionMessage(e)
        ),
        class = "ucminf_error"
      )
    }
  )
  
  out
}