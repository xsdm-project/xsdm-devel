# ============================================================================
# Graphical (not unit) tests for interpret_parameters() and auto_plot_lims_()
#
# Run interactively and verify each plot visually. To save plots to a PDF
# for later review, set save_pdf <- TRUE and specify a filename.
# ============================================================================

save_pdf <- FALSE
if (save_pdf) pdf("interpret_parameters_test_plots.pdf", width = 10, height = 6)

library(xsdm)

# ----------------------------------------------------------------------------
# Helper: safely call internal function
# ----------------------------------------------------------------------------
auto_plot_lims_ <- xsdm:::auto_plot_lims_

# ----------------------------------------------------------------------------
# 1) Legacy single-panel calls (no env_dat / occ)
# ----------------------------------------------------------------------------

param_list <- list(
  mu      = c(1, 2),
  sigltil = c(1, 3),
  sigrtil = c(2, 1),
  ctil    = 1,
  pd      = 0.5,
  o_mat   = diag(2)
)

# Two univariate plots in each direction
interpret_parameters(param_list, plot_indices = 1, plot_lims = list(c(-2, 4)))
interpret_parameters(param_list, plot_indices = 2, plot_lims = list(c(-6, 2)))

# Bivariate plot
interpret_parameters(
  param_list,
  plot_indices = c(1, 2),
  plot_lims    = list(c(-8, 16), c(-16, 8))
)

# Doubling limits: same contour shape, different axis ticks
interpret_parameters(
  param_list,
  plot_indices = c(1, 2),
  plot_lims    = list(c(-16, 32), c(-32, 16))
)

# Rotated parameters
theta <- 30 * pi / 180
param_list$o_mat <- matrix(
  c(cos(theta), sin(theta), -sin(theta), cos(theta)),
  nrow = 2, ncol = 2
)
interpret_parameters(param_list, plot_indices = 1, plot_lims = list(c(-2, 4)))
interpret_parameters(param_list, plot_indices = 2, plot_lims = list(c(-6, 2)))
interpret_parameters(
  param_list,
  plot_indices = c(1, 2),
  plot_lims    = list(c(-8, 16), c(-16, 8))
)

# ----------------------------------------------------------------------------
# 2) New two-panel behaviour (using package example data)
# ----------------------------------------------------------------------------

# Load example data (adjust names as needed; these are placeholders)
env_dat <- xsdm::example_1_env_array
occ     <- xsdm::example_1_occurrence_vector
pl      <- xsdm::example_1_param_list_example

# 2a) Bivariate, default breadth = 1 (full range + 10% margin)
interpret_parameters(
  pl,
  plot_indices = c(1, 2),
  env_dat      = env_dat,
  occ          = occ
)

# 2b) Bivariate, breadth = 0.5 (halfway)
interpret_parameters(
  pl,
  plot_indices = c(1, 2),
  env_dat      = env_dat,
  occ          = occ,
  breadth      = 0.5
)

# 2c) Bivariate, breadth = 0 (pinprick around mu)
interpret_parameters(
  pl,
  plot_indices = c(1, 2),
  env_dat      = env_dat,
  occ          = occ,
  breadth      = 0
)

# 2d) Bivariate with user-supplied plot_lims (using internal helper)
lims_user <- auto_plot_lims_(env_dat, pl, c(1, 2), breadth = 1)
interpret_parameters(
  pl,
  plot_indices = c(1, 2),
  plot_lims    = lims_user,
  env_dat      = env_dat,
  occ          = occ
)

# 2e) Univariate two-panel
interpret_parameters(pl, plot_indices = 1, env_dat = env_dat, occ = occ)
interpret_parameters(pl, plot_indices = 2, env_dat = env_dat, occ = occ)

# 2f) Univariate narrower breadth
interpret_parameters(
  pl,
  plot_indices = 1,
  env_dat      = env_dat,
  occ          = occ,
  breadth      = 0.3
)

# ----------------------------------------------------------------------------
# 3) auto_plot_lims_ standalone tests
# ----------------------------------------------------------------------------

# 3a) breadth = 1: limits cover all env values
lims_full <- auto_plot_lims_(env_dat, pl, c(1, 2), breadth = 1)
for (k in 1:2) {
  env_k <- as.numeric(env_dat[, , k])
  lo <- pl$mu[k] + lims_full[[k]][1]
  hi <- pl$mu[k] + lims_full[[k]][2]
  stopifnot(lo < min(env_k))
  stopifnot(hi > max(env_k))
  message(sprintf(
    "Var %d: data [%.2f, %.2f] inside plot window [%.2f, %.2f] -- OK",
    k, min(env_k), max(env_k), lo, hi
  ))
}

# 3b) breadth = 0: pinprick width < 1e-4
lims_pin <- auto_plot_lims_(env_dat, pl, c(1, 2), breadth = 0)
for (k in 1:2) {
  width <- lims_pin[[k]][2] - lims_pin[[k]][1]
  stopifnot(width < 1e-4)
  message(sprintf("Var %d: pinprick width %.2e -- OK", k, width))
}

# 3c) Monotonicity: larger breadth gives wider limits
lims_lo <- auto_plot_lims_(env_dat, pl, c(1, 2), breadth = 0.2)
lims_hi <- auto_plot_lims_(env_dat, pl, c(1, 2), breadth = 0.8)
for (k in 1:2) {
  w_lo <- lims_lo[[k]][2] - lims_lo[[k]][1]
  w_hi <- lims_hi[[k]][2] - lims_hi[[k]][1]
  stopifnot(w_hi >= w_lo)
  message(sprintf("Var %d: width 0.2 = %.4f < width 0.8 = %.4f -- OK",
                  k, w_lo, w_hi))
}

# 3d) Test margin parameter (default 0.1, try 0.2)
lims_margin <- auto_plot_lims_(env_dat, pl, c(1, 2), breadth = 1, margin = 0.2)
for (k in 1:2) {
  env_range <- diff(range(env_dat[, , k]))
  default_margin <- 0.1 * env_range
  new_margin <- 0.2 * env_range
  # The width (hi - lo) should be increased by extra 0.2 * range compared to default?
  # Actually the extra margin is added symmetrically, so total width increases by 2 * extra_margin.
  # But for simplicity, just check that limits are wider:
  lims_default <- auto_plot_lims_(env_dat, pl, c(1, 2), breadth = 1, margin = 0.1)
  stopifnot(lims_margin[[k]][1] < lims_default[[k]][1])
  stopifnot(lims_margin[[k]][2] > lims_default[[k]][2])
  message(sprintf("Var %d: margin 0.2 gives wider limits than margin 0.1 -- OK", k))
}

# ----------------------------------------------------------------------------
# 4) Expected-error cases
# ----------------------------------------------------------------------------

# 4a) Neither plot_lims nor env_dat supplied
tryCatch(
  interpret_parameters(param_list, plot_indices = 1),
  error = function(e) message("OK, got expected error:\n  ", conditionMessage(e))
)

# 4b) occ without env_dat
tryCatch(
  interpret_parameters(param_list, plot_indices = 1, occ = occ),
  error = function(e) message("OK, got expected error:\n  ", conditionMessage(e))
)

# 4c) breadth out of range
tryCatch(
  interpret_parameters(pl, c(1, 2), env_dat = env_dat, occ = occ, breadth = 1.5),
  error = function(e) message("OK, got expected error:\n  ", conditionMessage(e))
)
tryCatch(
  interpret_parameters(pl, c(1, 2), env_dat = env_dat, occ = occ, breadth = -0.01),
  error = function(e) message("OK, got expected error:\n  ", conditionMessage(e))
)

# 4d) Invalid plot_indices (out of range)
tryCatch(
  interpret_parameters(pl, plot_indices = 3, env_dat = env_dat, occ = occ),
  error = function(e) message("OK, got expected error:\n  ", conditionMessage(e))
)

# 4e) Non-finite mu
bad_param <- pl
bad_param$mu[1] <- NA
tryCatch(
  interpret_parameters(bad_param, plot_indices = 1, plot_lims = list(c(-2,4))),
  error = function(e) message("OK, got expected error:\n  ", conditionMessage(e))
)

# ----------------------------------------------------------------------------
# 5) Additional graphical argument passing test
# ----------------------------------------------------------------------------
# Check that ... is passed correctly to plot/image
interpret_parameters(
  pl,
  plot_indices = 1,
  env_dat = env_dat,
  occ = occ,
  breadth = 0.5,
  col.main = "blue",    # passed to title()
  cex.lab = 1.2
)

# ----------------------------------------------------------------------------

if (save_pdf) dev.off()
message("\nAll tests completed. Please review the generated plots manually.")