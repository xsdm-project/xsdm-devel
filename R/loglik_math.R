#' Log-likelihood function for the xsdm model, parameters on the math scale.
#'
#' Computes the log-likelihood for the xsdm model given environmental data, a
#' vector of occurrences and pseudo-absences, and model parameters on the
#' math scale. This is the function that one optimizes to fit xsdm with data.
#'
#' @param param_vector A **named numeric vector** of math-scale parameters.
#' When \code{mask = NULL}, the names must exactly match the canonical schema
#' returned by \code{\link{make_mask_names}(p)} where \code{p = dim(env_dat)[3]},
#' and the length must equal \code{\link{num_par}(p)}.
#' When \code{mask} is supplied, \code{param_vector} should contain only the
#' names **not** present in \code{mask}, in the canonical order.
#' In both cases the vector is combined with \code{mask} via
#' \code{\link{create_param_vector_masked}} and then mapped to biological-scale
#' parameters via \code{\link{math_to_bio}}.  Must not contain missing values.
#' See Details for the full naming and ordering conventions.
#' @param env_dat The environmental data array, dimensions
#' \code{n_loc x n_time x p} (number of locations x time-series length x number
#' of environmental variables). Must be a 3-dimensional array with no missing
#' values.
#' @param occ Presence/pseudo-absence binary vector. Same length as dimension 1 of
#' \code{env_dat}.
#' @param mask For optionally keeping some parameters at fixed values
#' during optimization. Either NULL or a named numeric vector. The NULL case
#' means all parameters are in param_vector, corresponding to the case where all
#' parameters will be adjustable by the optimizer when this function is passed
#' to it as the objective function. In the non-NULL case, names of entries of
#' \code{mask} must correspond precisely to parameter names (see Details),
#' and then those values are used. In that case, \code{param_vector} is construed
#' to contain the values of the other parameters, in order (see Details). So, in
#' particular, the length of \code{param_vector} plus the length of
#' \code{mask} must equal the total number of parameters for the model, which is
#' determined by dim(env_dat)[3]. The most common case for most applications
#' will be \code{mask=NULL}. Entries of \code{mask} are interpreted on the math
#' scale.
#' @param negative Logical. If TRUE returns the negative of the log-likelihood
#' instead of the log-likelihood itself. Facilitates optimization with some
#' optimizers.
#' @param num_threads Number of threads for parallel computation. Defaults to
#' \code{RcppParallel::defaultNumThreads()}.
#'
#' @returns A single value, the log-likelihood (or the negative log-likelihood,
#' if \code{negative} is TRUE).
#' @export
#'
#' @section Implementation:
#' This is a thin R wrapper around the C++ implementation
#' \code{loglik_math_cpp}; the optimizer hot path is pure C++. A pure-R
#' reference, \code{loglik_math_r}, is kept internal to the package and is
#' used only by the parity tests in
#' \code{tests/testthat/test-loglik_math_r_vs_cpp.R}.
#'
#' @details
#' Optimizing the likelihood and profiling requires conventions for transforming
#' parameters from unconstrained spaces to the constrained space of possible
#' parameters which can be accepted by \code{loglik_bio}. This function and
#' \code{math_to_bio} implement those conventions, and also allow for
#' optimizations while keeping one or more parameters fixed, including
#' potentially at boundary values. Typically \code{loglik_math} is the function
#' one optimizes numerically
#' in order to fit xsdm or a boundary model with data, or to profile a fitted
#' model. For what follows, denote \code{dim(env_dat)[3]} by \code{p}.
#'
#' We start by explaining the case \code{mask=NULL}, for which all model
#' parameters are in \code{param_vector}. The parameters of \code{param_vector}
#' are assumed to appear in the following order:
#' \enumerate{
#' \item Parameters for \code{mu}, of which there are \code{p};
#' \item Parameters which are \code{exp}-transformed to get the entries of
#' \code{sigltil}, of which there are \code{p};
#' \item Parameters which are \code{exp}-transformed to get the entries of
#' \code{sigrtil}, of which there are \code{p};
#' \item The parameter \code{ctil};
#' \item A parameter which is \code{expit} transformed to get \code{pd};
#' \item Parameters which are inserted via column-major order into the lower-
#' triangle of a skew-symmetric matrix which is then transformed by the matrix
#' exponential to get \code{o_mat}, of which there are (p^2-p)/2.
#' }
#' Thus, when \code{mask} is \code{NULL}, \code{param_vector} must be an
#' unconstrained numeric vector of length \code{3*p+2+(p^2-p)/2} with no missing
#' values.
#'
#' The argument \code{mask} is used in the event one wants to fix certain
#' parameters and optimize over the remaining parameters. This argument must be
#' a named numeric vector with unique names being some but not all of the
#' \code{3*p+2+(p^2-p)/2} following: \code{mu1}, \code{mu2}, \ldots, \code{mup},
#' \code{sigltil1}, \code{sigltil2}, \ldots, \code{sigltilp}, \code{sigrtil1},
#' \code{sigrtil2}, \ldots, \code{sigrtilp}, \code{ctil}, \code{pd}, and
#' \code{o_mati} for \code{i} ranging from 1 to \code{(p^2-p)/2}. These names
#' must be used exactly. See the function \code{make_mask_names}, which
#' facilitates the construction of a correctly formatted \code{mask} argument.
#' Entries of \code{mask} are on the math scale; \code{bio_to_math} can convert
#' biological-scale constraints to math scale.
#'
#' The missing entries of \code{mask} are filled in using the entries of
#' \code{param_vector}, in the order specified above, and then the
#' transformations described above (implemented by \code{math_to_bio}) are applied
#' to get biological-scale parameters which are passed to \code{loglik_bio} to
#' get the log likelihood.
#'
#' Entries of \code{mask} corresponding to \code{sigltil} or \code{sigrtil}
#' can be \code{Inf}.  Likewise, the entry of \code{mask} corresponding
#' to \code{pd} can be \code{Inf} (on the math scale, corresponding to a
#' biological-scale value of 1). This functionality is used to fit boundary
#' models. Entries of \code{param_vector} must be finite.
#'
#' @examples
#' # Testing the function with the example data
#' loglik_math(
#'   param_vector = example_1$par_vec,
#'   env_dat = example_1$env_array,
#'   occ = example_1$occ_vec
#' )
#' # Mute one parameter to use the mask
#' par_vec <- example_1$par_vec[-2]
#' mask_parameters_a <- c(mu2 = 6.5)
#' loglik_math(
#'   param_vector = par_vec,
#'   env_dat = example_1$env_array,
#'   occ = example_1$occ_vec,
#'   mask = mask_parameters_a
#' )
#' # Return the negative
#' loglik_math(
#'   param_vector = example_1$par_vec,
#'   env_dat = example_1$env_array,
#'   occ = example_1$occ_vec,
#'   negative = TRUE
#' )

loglik_math <- function(param_vector,
                        env_dat,
                        occ,
                        mask = NULL,
                        num_threads = RcppParallel::defaultNumThreads(),
                        negative = TRUE) {
  # Input validation (kept in R so the user-facing error messages and
  # checkmate-style asserts continue to match the historical API surface;
  # the math itself runs entirely in C++ via `loglik_math_cpp`).

  checkmate::assert(
    checkmate::check_logical(occ, any.missing = FALSE),
    checkmate::check_integerish(occ, lower = 0, upper = 1, any.missing = FALSE),
    .var.name = "occ"
  )

  check_env_array(env_dat)

  param_vector <- unlist(param_vector)
  checkmate::assert_vector(param_vector, any.missing = FALSE)
  p <- dim(env_dat)[3]

  # When `param_vector` arrives unnamed (e.g. from a C++ optimizer callback),
  # assign canonical free-parameter names so the C++ side can validate and
  # overlay them correctly.
  if (is.null(names(param_vector))) {
    all_canonical <- names(make_mask_names(p))
    free_names    <- if (!is.null(mask)) setdiff(all_canonical, names(mask)) else all_canonical
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

  loglik_math_cpp(
    param_vector = param_vector,
    env_dat      = env_dat,
    occ          = as.integer(occ),
    mask         = mask,
    negative     = isTRUE(negative),
    num_threads  = as.integer(num_threads)
  )
}
