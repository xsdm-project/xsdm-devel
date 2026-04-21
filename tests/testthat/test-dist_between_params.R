library(testthat)
test_that("dist_between_params: basic structure, correctness,
          and parity with slow version", {
  skip_if_not_installed("clue")
  skip_if_not_installed("expm")
  
  # ---- Simple 2D deterministic setup -----------------------------------------
  # Rotation matrices near each other, same sigs and other scalars
  ang1 <- -46 * pi / 180
  o_mat1 <- matrix(c(cos(ang1), sin(ang1), -sin(ang1), cos(ang1)), 2, 2)
  sigltil1 <- c(1, 2)
  sigrtil1 <- c(3, 4)
  mu <- c(0, 0)
  ctil <- 10
  pd <- 0.98
  p1 <- list(
    mu = mu,
    ctil = ctil,
    pd = pd,
    o_mat = o_mat1,
    sigltil = sigltil1,
    sigrtil = sigrtil1
  )
  
  ang2 <- -44 * pi / 180
  o_mat2 <- matrix(c(cos(ang2), sin(ang2), -sin(ang2), cos(ang2)), 2, 2)
  sigltil2 <- c(1, 2)
  sigrtil2 <- c(3, 4)
  p2 <- list(
    mu = mu,
    ctil = ctil,
    pd = pd,
    o_mat = o_mat2,
    sigltil = sigltil2,
    sigrtil = sigrtil2
  )
  
  # ---- Return types and names -------------------------------------------------
  res1 <- dist_between_params(p1, p2)
  res2 <- dist_between_params(p1, p2, give_closest_rep = TRUE)
  
  expect_type(res1, "double")
  expect_length(res1, 1)
  expect_type(res2, "list")
  expect_named(res2, c("distance", "representative"))
  expect_equal(res1, res2$distance)
  
  rep <- res2$representative
  expect_type(rep, "list")
  expect_named(rep, c("mu", "sigltil", "sigrtil", "ctil", "pd", "o_mat"))
  expect_equal(rep$mu, mu)
  expect_equal(rep$ctil, ctil)
  expect_equal(rep$pd, pd)
  
  # ---- Accuracy for this setup ------------------------------------------------
  # Here sig terms and scalars match exactly, so distance reduces to Frobenius
  # distance between the o_mat columns under optimal assignment/sign.
  expect_equal(res1, sqrt(sum((o_mat1 - o_mat2)^2)))
  
  # ---- Effect of changing the non-equivalence-class parameters ----------------
  p2b <- p2
  p2b$mu <- c(1, 1)
  p2b$ctil <- 3
  p2b$pd <- 0.1
  res2b <- dist_between_params(p1, p2b, give_closest_rep = TRUE)
  expect_equal(
    res2b$distance,
    sqrt(res1^2 + (1^2 + 1^2) + (10 - 3)^2 + (0.98 - 0.1)^2)
  )
  
  # ---- Parity with the slow exhaustive version (3D) ---------------------------
  set.seed(123)
  S3 <- matrix(rnorm(9), 3, 3)
  S3 <- S3 - t(S3)
  o_mat <- expm::expm(S3)
  sigltil <- runif(3)
  sigrtil <- runif(3)
  mu <- runif(3)
  ctil <- rnorm(1)
  pd <- runif(1)
  p1 <- list(mu = mu, sigltil = sigltil, sigrtil = sigrtil, ctil = ctil, pd = pd, o_mat = o_mat)
  
  S3 <- matrix(rnorm(9), 3, 3)
  S3 <- S3 - t(S3)
  o_mat <- expm::expm(S3)
  sigltil <- runif(3)
  sigrtil <- runif(3)
  mu <- runif(3)
  ctil <- rnorm(1)
  pd <- runif(1)
  p2 <- list(mu = mu, sigltil = sigltil, sigrtil = sigrtil, ctil = ctil, pd = pd, o_mat = o_mat)
  
  res_slow <- xsdm:::distance_between_params(p1, p2, GiveClosestRep = TRUE)
  res_fast <- dist_between_params(p1, p2, give_closest_rep = TRUE)
  expect_equal(res_slow$distance, res_fast$distance, tolerance = 1e-12)
  expect_equal(res_slow$representative$mu, res_fast$representative$mu)
  expect_equal(res_slow$representative$ctil, res_fast$representative$ctil)
  expect_equal(res_slow$representative$pd, res_fast$representative$pd)
  expect_equal(res_slow$representative$sigltil, res_fast$representative$sigltil)
  expect_equal(res_slow$representative$sigrtil, res_fast$representative$sigrtil)
  expect_equal(res_slow$representative$o_mat, res_fast$representative$o_mat, tolerance = 1e-12)
  
  # ---- Parity with the slow exhaustive version (4D) ---------------------------
  set.seed(456)
  S4 <- matrix(rnorm(16), 4, 4)
  S4 <- S4 - t(S4)
  o_mat <- expm::expm(S4)
  sigltil <- runif(4)
  sigrtil <- runif(4)
  mu <- runif(4)
  ctil <- rnorm(1)
  pd <- runif(1)
  p1 <- list(mu = mu, sigltil = sigltil, sigrtil = sigrtil, ctil = ctil, pd = pd, o_mat = o_mat)
  
  S4 <- matrix(rnorm(16), 4, 4)
  S4 <- S4 - t(S4)
  o_mat <- expm::expm(S4)
  sigltil <- runif(4)
  sigrtil <- runif(4)
  mu <- runif(4)
  ctil <- rnorm(1)
  pd <- runif(1)
  p2 <- list(mu = mu, sigltil = sigltil, sigrtil = sigrtil, ctil = ctil, pd = pd, o_mat = o_mat)
  
  res_slow <- xsdm:::distance_between_params(p1, p2, GiveClosestRep = TRUE)
  res_fast <- dist_between_params(p1, p2, give_closest_rep = TRUE)
  expect_equal(res_slow$distance, res_fast$distance, tolerance = 1e-12)
  expect_equal(res_slow$representative$mu, res_fast$representative$mu)
  expect_equal(res_slow$representative$sigltil, res_fast$representative$sigltil)
  expect_equal(res_slow$representative$sigrtil, res_fast$representative$sigrtil)
  expect_equal(res_slow$representative$ctil, res_fast$representative$ctil)
  expect_equal(res_slow$representative$pd, res_fast$representative$pd)
  expect_equal(res_slow$representative$o_mat, res_fast$representative$o_mat, tolerance = 1e-12)
  
  # ---- Equivalence of math vs biological input formats -----------------------
  # p = 1, 2, 3, 4 cases
  for (pp in 1:4) {
    set.seed(101 + pp)
    nm <- make_mask_names(pp)
    v1 <- nm
    v2 <- nm
    v1[seq_along(v1)] <- rnorm(length(v1))
    v2[seq_along(v2)] <- rnorm(length(v2))
    
    p1_bio <- math_to_bio(v1)
    p2_bio <- math_to_bio(v2)
    
    d_math <- dist_between_params(v1, v2)
    d_bio <- dist_between_params(p1_bio, p2_bio)
    expect_equal(d_math, d_bio, tolerance = 1e-10)
  }
  
  # ---- Math-scale stress case: very large sigma value -------------------------
  set.seed(202)
  nm <- make_mask_names(1)
  v1 <- nm
  v2 <- nm
  v1[seq_along(v1)] <- rnorm(length(v1))
  v2[seq_along(v2)] <- rnorm(length(v2))
  v1["sigltil1"] <- 1000
  v2["sigltil1"] <- 1001
  expect_true(is.finite(dist_between_params(v1, v2)))
  
  # ---- Mask complementarity case (math-scale, with Inf allowed only in sig*/pd)
  nm <- make_mask_names(2)
  set.seed(303)
  v <- nm
  v[seq_along(v)] <- rnorm(length(v))
  # Remove one entry from the free vector, put it (as Inf) into the mask
  # Here we treat sigltil1 as masked/Inf; the union (v + mask) still covers canon.
  mask <- c(sigltil1 = Inf)
  v_drop <- setdiff(names(v), names(mask))[1]
  v_masked <- v[names(v) != names(mask)]
  # symmetrical case to itself should be zero distance
  expect_equal(dist_between_params(p1 = v_masked,
                                   p2 = v_masked,
                                   mask = mask),
               0)
})


test_that("dist_between_params: core behavior, types, and basic accuracy", {
  # 2D deterministic setup: small angle difference, same sigs/pd/ctil
  ang1 <- -46 * pi / 180
  o_mat1 <- matrix(c(cos(ang1), sin(ang1), -sin(ang1), cos(ang1)), 2, 2)
  sigltil1 <- c(1, 2)
  sigrtil1 <- c(3, 4)
  mu <- c(0, 0); ctil <- 10; pd <- 0.98
  p1 <- list(mu = mu, ctil = ctil, pd = pd, o_mat = o_mat1,
             sigltil = sigltil1, sigrtil = sigrtil1)
  
  ang2 <- -44 * pi / 180
  o_mat2 <- matrix(c(cos(ang2), sin(ang2), -sin(ang2), cos(ang2)), 2, 2)
  p2 <- list(mu = mu, ctil = ctil, pd = pd, o_mat = o_mat2,
             sigltil = sigltil1, sigrtil = sigrtil1)
  
  # Return types
  res1 <- dist_between_params(p1, p2)
  res2 <- dist_between_params(p1, p2, give_closest_rep = TRUE)
  expect_type(res1, "double"); expect_length(res1, 1)
  expect_type(res2, "list"); expect_named(res2, c("distance", "representative"))
  expect_equal(res1, res2$distance)
  
  # Representative has canonical names and copied scalars
  rep <- res2$representative
  expect_named(rep, c("mu","sigltil","sigrtil","ctil","pd","o_mat"))
  expect_equal(rep$mu, mu)
  expect_equal(rep$ctil, ctil); expect_equal(rep$pd, pd)
  
  # With identical sigs/scalars, distance reduces to Frobenius distance on o_mat
  # after best assignment/sign
  expect_equal(res1, sqrt(sum((o_mat1 - o_mat2)^2)))
})

test_that("dist_between_params: tie branch (pos == neg) sets sign to +1", {
  # Construct a tie for one column: make o_mat columns orthogonal (a·b = 0) and
  # set sigltil==sigrtil in both p1 and p2 so sig terms match for pos/neg.
  # Then ||a-b||^2 == ||a+b||^2, yielding a tie; code should pick +1.
  
  # p = 2
  o1 <- diag(2)            # columns: e1, e2
  o2 <- matrix(c(0,1,1,0), 2, 2)  # columns: e2, e1 (orthogonal swap)
  mu <- c(0,0); ctil <- 0.5; pd <- 0.3
  sig_equal <- c(2, 2)     # sigltil == sigrtil for tie in sig terms
  
  p1 <- list(mu = mu, ctil = ctil, pd = pd, o_mat = o1,
             sigltil = sig_equal, sigrtil = sig_equal)
  p2 <- list(mu = mu, ctil = ctil, pd = pd, o_mat = o2,
             sigltil = sig_equal, sigrtil = sig_equal)
  
  # We can't directly access posneg, but we can check that the chosen
  # representative corresponds to +1 flips for the paired columns
  # (tie -> +1 per code path).
  out <- dist_between_params(p1, p2, give_closest_rep = TRUE)
  rep <- out$representative
  
  # In this setup, the optimal assignment pairs columns {1<->2, 2<->1} with + sign.
  # So rep$o_mat should equal p1$o_mat with columns swapped but not sign-flipped.
  expect_equal(rep$o_mat, p1$o_mat[, c(2,1), drop = FALSE])
  # sigltil/sigrtil should also be swapped accordingly, not interchanged between sigl/sigr
  expect_equal(rep$sigltil, p1$sigltil[c(2,1)])
  expect_equal(rep$sigrtil, p1$sigrtil[c(2,1)])
})


test_that("dist_between_params: math-scale inputs + mask complementarity
          (including Inf for sig*/pd)", {
  # p = 2; build full canonical name vector, then split between p and mask
  nm <- make_mask_names(2)
  set.seed(303)
  # full canonical math-scale vector (length 9)
  v <- nm; v[seq_along(v)] <- rnorm(length(v))  
  
  # ---- Valid mask branch (control) ----
  # allowed names for Inf
  mask <- c(sigltil1 = Inf, pd = Inf)                 
  # drop the masked names from the free vector
  v_free <- v[setdiff(names(v), names(mask))]         
  d0 <- dist_between_params(p1 = v_free, p2 = v_free, mask = mask)
  expect_equal(d0, 0)
  
  # ---- Invalid mask branch (what we want to test) ----
  # To reach the "mask can only have infinite values ..." check,
  # ensure the union of names matches the canonical set: remove the masked name
  # from the free vector and put it in 'mask'.
  bad_mask <- c(mu1 = Inf)                            # NOT allowed to be Inf
  # remove mu1 so |v_free_bad| + |mask| == 9
  v_free_bad <- v[setdiff(names(v), names(bad_mask))] 
  expect_error(
    dist_between_params(p1 = v_free_bad, p2 = v_free_bad, mask = bad_mask),
    regexp = "Must be TRUE"
  )
})


test_that("dist_between_params:
          mismatched lengths/names on math-scale raise errors", {
  # --- 1) Wrong length: assert_integerish() ---
  nm <- make_mask_names(2)
  v_bad_len <- nm[-1]; v_bad_len[] <- 0
  expect_error(
    dist_between_params(v_bad_len, v_bad_len),
    regexp = "Must be of type 'integerish'",
    fixed  = FALSE
  )
  
  # --- 2) Wrong union of names: assert_true() => "Must be TRUE" ---
  nm <- make_mask_names(2)
  v  <- nm; v[] <- 0
  
  # Remove *mu2* from the free vector (length 8)
  v_free_bad_union <- v[names(v) != "mu2"]
  
  # Put *mu1* in the mask (length 1) -> total remains 9, but the union of names
  # is canonical minus "mu2" plus a duplicate "mu1", so setequal(...) is FALSE.
  mask_overlap <- c(mu1 = 0.1)
  
  expect_error(
    dist_between_params(v_free_bad_union,
                        v_free_bad_union,
                        mask = mask_overlap),
    regexp = "Must be TRUE",
    fixed  = FALSE
  )
})



test_that("dist_between_params", {
  skip_if_not_installed("clue")
  skip_if_not_installed("expm")
  
  set.seed(123)
  S3 <- matrix(rnorm(9), 3, 3); S3 <- S3 - t(S3)
  o1 <- expm::expm(S3)
  p1 <- list(mu = runif(3), sigltil = runif(3), sigrtil = runif(3),
             ctil = rnorm(1), pd = runif(1), o_mat = o1)
  
  S3 <- matrix(rnorm(9), 3, 3); S3 <- S3 - t(S3)
  o2 <- expm::expm(S3)
  p2 <- list(mu = runif(3), sigltil = runif(3), sigrtil = runif(3),
             ctil = rnorm(1), pd = runif(1), o_mat = o2)
  
  res_slow <- xsdm:::distance_between_params(p1, p2, GiveClosestRep = TRUE)
  res_fast <- dist_between_params(p1, p2, give_closest_rep = TRUE)
  expect_equal(res_slow$distance, res_fast$distance, tolerance = 1e-12)
  expect_equal(res_slow$representative$o_mat, res_fast$representative$o_mat,
               tolerance = 1e-12)
})
