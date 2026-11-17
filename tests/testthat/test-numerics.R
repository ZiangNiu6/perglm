test_that("cached residualization agrees with weighted least squares", {
  z <- cbind(
    intercept = 1,
    z1 = c(-2, -1, 0, 1, 2, 3),
    z2 = c(0, 1, 0, -1, 1, 2)
  )
  weights <- c(0.2, 1, 4, 0.5, 2, 3)
  projection <- prepare_projection(z, weights)

  for (exposure in list(c(0, 1, 0, 0, 1, 0), c(2, -3, 1, 5, 0, 4))) {
    residuals <- residualize_exposure(projection, exposure)
    reference <- stats::lm.wfit(z, exposure, weights)$residuals
    expect_equal(residuals, reference, tolerance = 1e-12)
    expect_equal(
      unname(drop(crossprod(z, weights * residuals))),
      rep(0, ncol(z)),
      tolerance = 1e-11
    )
  }
})

test_that("projection respects nuisance adjustment and common weight scaling", {
  z <- cbind(1, c(-2, -1, 0, 1, 2, 3))
  exposure <- c(1, 0, 0, 2, 0, 3)
  weights <- c(1, 2, 3, 0.5, 4, 1)
  residuals <- residualize_exposure(prepare_projection(z, weights), exposure)

  for (scale in c(1e-250, 1e250)) {
    expect_equal(
      residualize_exposure(prepare_projection(z, weights * scale), exposure),
      residuals,
      tolerance = 1e-12
    )
  }
  adjusted <- 2 * exposure + drop(z %*% c(3, -4))
  expect_equal(
    residualize_exposure(prepare_projection(z, weights), adjusted),
    2 * residuals,
    tolerance = 1e-12
  )
})

test_that("projection is stable with a wide range of positive weights", {
  z <- matrix(1, nrow = 4, ncol = 1)
  exposure <- c(-2, 1, 4, 7)
  weights <- c(1e-300, 1e-100, 1e100, 1e300)
  expect_equal(
    residualize_exposure(prepare_projection(z, weights), exposure),
    exposure - 7,
    tolerance = 1e-12
  )
})

test_that("invalid projection designs and tolerances are rejected", {
  z <- cbind(1, c(-1, 0, 1))
  expect_error(prepare_projection(cbind(1, 1, c(-1, 0, 1)), rep(1, 3)),
               "full column rank")
  expect_error(prepare_projection(matrix(1, 2, 3), rep(1, 2)),
               "at least as many rows")
  expect_error(prepare_projection(matrix(numeric(), 0, 1), numeric()),
               "nonempty")
  expect_error(prepare_projection(z * NA_real_, rep(1, 3)), "finite real")
  expect_error(prepare_projection(z + 0i, rep(1, 3)), "finite real")
  expect_error(prepare_projection(as.data.frame(z), rep(1, 3)), "matrix")
  for (tol in list(0, 1, NA_real_, Inf, c(0.1, 0.2), 1i)) {
    expect_error(prepare_projection(z, rep(1, 3), tol), "tol must")
  }
})

test_that("invalid weights and exposures are rejected", {
  z <- cbind(1, c(-1, 0, 1))
  for (weights in list(c(1, 1), c(1, 0, 1), c(1, -1, 1),
                       c(1, NA, 1), c(1, Inf, 1), rep(1i, 3),
                       matrix(1, 3, 1))) {
    expect_error(prepare_projection(z, weights), "weights must")
  }
  projection <- prepare_projection(z, rep(1, 3))
  for (exposure in list(c(0, 1), c(0, NA, 1), c(0, Inf, 1),
                        rep(1i, 3), matrix(1, 3, 1))) {
    expect_error(residualize_exposure(projection, exposure), "exposure must")
  }
  expect_error(residualize_exposure(NULL, c(0, 1, 0)), "projection must")
})

test_that("studentization uses the empirical variance with divisor n", {
  expect_equal(studentize_score(c(1, 2, 3)), 6 / sqrt(2), tolerance = 1e-14)
  expect_equal(studentize_score(c(-2, -1, 1, 2)), 0)
})

test_that("studentization is stable under extreme common scaling", {
  contributions <- c(-2, 1, 3, 5)
  reference <- sum(contributions) /
    sqrt(sum((contributions - mean(contributions))^2))
  for (scale in c(1e-250, 1, 1e250)) {
    expect_equal(studentize_score(scale * contributions), reference,
                 tolerance = 1e-14)
  }
  expect_equal(studentize_score(-contributions), -reference, tolerance = 1e-14)
})

test_that("studentization preserves tiny differences between large or tiny scores", {
  smallest <- .Machine$double.xmin * .Machine$double.eps
  for (contributions in list(c(1e308, 1e308 + 2e292),
                            c(1e-308, 1e-308 + smallest))) {
    scale <- max(abs(contributions))
    # For two scores, the centered norm is abs(r[2] - r[1]) / sqrt(2).
    expected <- sqrt(2) * sum(contributions / scale) /
      (abs(contributions[2] - contributions[1]) / scale)
    expect_equal(studentize_score(contributions), expected, tolerance = 1e-14)
  }
})

test_that("studentization handles anchor differences that overflow", {
  contributions <- c(1e308, -1e308, 1e308)
  expect_equal(studentize_score(contributions), sqrt(3 / 8), tolerance = 1e-14)
})

test_that("studentization rejects invalid and degenerate contributions", {
  for (contributions in list(numeric(), 1, c(1, NA), c(1, Inf),
                            c(1, NaN), c(1i, 2i), matrix(1:4, 2, 2))) {
    expect_error(studentize_score(contributions), "contributions must")
  }
  for (contributions in list(c(0, 0), rep(2, 4), rep(-2, 4))) {
    expect_error(studentize_score(contributions), "positive empirical variance")
  }
})
