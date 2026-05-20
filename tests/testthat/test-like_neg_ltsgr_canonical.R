library(testthat)
library(xsdm)

# ---------------------------------------------------------------------------
# Confirms the canonical rename: like_neg_ltsgr is the exported public name
# and like_neg_ltsgr_cpp is the unexported back-compat alias. Both must
# return identical values on the same inputs.
# ---------------------------------------------------------------------------

test_that("like_neg_ltsgr is exported and like_neg_ltsgr_cpp is unexported", {
  expect_true(exists("like_neg_ltsgr",
                     envir = asNamespace("xsdm"), inherits = FALSE))
  expect_false("like_neg_ltsgr_cpp" %in% getNamespaceExports("xsdm"))
  expect_true("like_neg_ltsgr"      %in% getNamespaceExports("xsdm"))
})

test_that("canonical and alias return identical values", {
  out_canonical <- like_neg_ltsgr(
    env_dat = example_1$env_array,
    mu      = example_1$true_par_list$mu,
    sigltil = example_1$true_par_list$sigltil,
    sigrtil = example_1$true_par_list$sigrtil,
    o_mat   = example_1$true_par_list$o_mat,
    num_threads = 1L
  )
  out_alias <- xsdm:::like_neg_ltsgr_cpp(
    env_dat = example_1$env_array,
    mu      = example_1$true_par_list$mu,
    sigltil = example_1$true_par_list$sigltil,
    sigrtil = example_1$true_par_list$sigrtil,
    o_mat   = example_1$true_par_list$o_mat,
    num_threads = 1L
  )
  expect_identical(out_canonical, out_alias)
})
