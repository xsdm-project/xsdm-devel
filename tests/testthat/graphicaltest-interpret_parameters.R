library(xsdmMle)
# Graphical (not unit) tests for interpret_parameters(). Run interactively
# and verify each plot visually.
#
# Section 1  Legacy single-panel (no env_dat / occ). Must be identical to
#            the pre-#44 output.
# Section 2  New two-panel (#44) using the package example data.
#            Axis limits auto-derived from the full env range + 10% margin.
# Section 3  auto_plot_lims stand-alone smoke tests.
# Section 4  Expected-error cases.

# ========================================================================
# 1) Legacy single-panel calls
# ========================================================================

param_list <- list(
  mu      = c(1, 2),
  sigltil = c(1, 3),
  sigrtil = c(2, 1),
  ctil    = 1,
  pd      = 0.5,
  o_mat   = diag(2)
)

# two univariate plots in each direction
interpret_parameters(param_list, plot_indices = 1, plot_lims = list(c(-2, 4)))
interpret_parameters(param_list, plot_indices = 2, plot_lims = list(c(-6, 2)))

# bivariate plot
interpret_parameters(
  param_list,
  plot_indices = c(1, 2),
  plot_lims    = list(c(-8, 16), c(-16, 8))
)

# doubling limits: same contour shape, different axis ticks
interpret_parameters(
  param_list,
  plot_indices = c(1, 2),
  plot_lims    = list(c(-16, 32), c(-32, 16))
)

# rotated parameters
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

# ========================================================================
# 2) New two-panel behaviour (#44)
#    Axis limits are derived from the FULL env range (all locations and
#    time steps, not filtered by occ) plus a 10% margin on each side.
#    The `breadth` knob controls the width with the same semantics as in
#    get_range_df():
#      breadth = 1  ->  full min-max env range + 10% margin
#      breadth = 0  ->  pinprick around mu
# ========================================================================

env_dat <- example_1_env_array
occ     <- example_1_occurrence_vector
pl      <- example_1_param_list_example

# 2a) Bivariate, default breadth = 1 (full range + margin).
#     EXPECTED: all scatter points (presences left, non-detections right)
#     sit comfortably inside the axis limits with breathing room on every
#     edge.
interpret_parameters(
  pl,
  plot_indices = c(1, 2),
  env_dat      = env_dat,
  occ          = occ
)

# 2b) Bivariate, breadth = 0.5 (halfway). Axis range visibly narrower than
#     2a; some data points near the tails may now sit close to the plot edge
#     or outside it.
interpret_parameters(
  pl,
  plot_indices = c(1, 2),
  env_dat      = env_dat,
  occ          = occ,
  breadth      = 0.5
)

# 2c) Bivariate, breadth = 0 (pinprick). Essentially no data visible;
#     smoke check that the code path does not blow up.
interpret_parameters(
  pl,
  plot_indices = c(1, 2),
  env_dat      = env_dat,
  occ          = occ,
  breadth      = 0
)

# 2d) Bivariate with user-supplied plot_lims. When plot_lims is given,
#     breadth is silently ignored; env_dat / occ still trigger two-panel.
interpret_parameters(
  pl,
  plot_indices = c(1, 2),
  plot_lims    = auto_plot_lims(env_dat, pl, c(1, 2)),
  env_dat      = env_dat,
  occ          = occ
)

# 2e) Univariate two-panel: growth curve + jitter cloud of env values at
#     presences (left) and non-detections (right). Same curve on both.
interpret_parameters(pl, plot_indices = 1, env_dat = env_dat, occ = occ)
interpret_parameters(pl, plot_indices = 2, env_dat = env_dat, occ = occ)

# 2f) Univariate narrower breadth.
interpret_parameters(
  pl,
  plot_indices = 1,
  env_dat      = env_dat,
  occ          = occ,
  breadth      = 0.3
)

# ========================================================================
# 3) auto_plot_lims standalone smoke tests
# ========================================================================

# 3a) breadth = 1: limits should comfortably cover all env values.
lims_full <- auto_plot_lims(env_dat, pl, c(1, 2), breadth = 1)
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

# 3b) breadth = 0: limits should collapse to a pinprick around mu.
lims_pin <- auto_plot_lims(env_dat, pl, c(1, 2), breadth = 0)
for (k in 1:2) {
  width <- lims_pin[[k]][2] - lims_pin[[k]][1]
  stopifnot(width < 1e-4)
  message(sprintf("Var %d: pinprick width %.2e -- OK", k, width))
}

# 3c) breadth monotonicity: wider breadth -> wider limits.
lims_lo <- auto_plot_lims(env_dat, pl, c(1, 2), breadth = 0.2)
lims_hi <- auto_plot_lims(env_dat, pl, c(1, 2), breadth = 0.8)
for (k in 1:2) {
  w_lo <- lims_lo[[k]][2] - lims_lo[[k]][1]
  w_hi <- lims_hi[[k]][2] - lims_hi[[k]][1]
  stopifnot(w_hi >= w_lo)
  message(sprintf("Var %d: width at 0.2 = %.4f < width at 0.8 = %.4f -- OK",
                  k, w_lo, w_hi))
}

# ========================================================================
# 4) Expected-error cases
# ========================================================================

# 4a) Neither plot_lims nor env_dat supplied.
tryCatch(
  interpret_parameters(param_list, plot_indices = 1),
  error = function(e) message("OK, got expected error:\n  ", conditionMessage(e))
)

# 4b) occ without env_dat.
tryCatch(
  interpret_parameters(param_list, plot_indices = 1, occ = occ),
  error = function(e) message("OK, got expected error:\n  ", conditionMessage(e))
)

# 4c) breadth out of range.
tryCatch(
  interpret_parameters(pl, c(1, 2), env_dat = env_dat, occ = occ, breadth = 1.5),
  error = function(e) message("OK, got expected error:\n  ", conditionMessage(e))
)
tryCatch(
  interpret_parameters(pl, c(1, 2), env_dat = env_dat, occ = occ, breadth = -0.01),
  error = function(e) message("OK, got expected error:\n  ", conditionMessage(e))
)
