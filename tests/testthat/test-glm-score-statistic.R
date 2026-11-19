.score_statistic_null <- function(family = gaussian()) {
  i <- seq_len(80)
  data <- data.frame(z = cos(i * 0.27), y = 0.3 + 0.2 * cos(i * 0.27) +
                       sin(i * 1.13) + 0.4 * cos(i * 0.83),
                     row.names = paste0("row", i))
  if (family$family == "poisson") {
    data$y <- (i * 7L %% 11L) + as.integer(data$z > 0)
  }
  fit_glm_null(y ~ z, data, family = family)
}

test_that("the observed GLM score agrees with independent R calculations", {
  for (family in list(gaussian(), poisson())) {
    null <- .score_statistic_null(family)
    for (exposure in list(sin(seq_len(null$nobs) * 0.39),
                          as.integer(seq_len(null$nobs) %% 3L == 0L))) {
      z <- null$projection$z
      w <- null$working_weights
      beta <- qr.coef(qr(sqrt(w) * z), sqrt(w) * exposure)
      r <- (exposure - drop(z %*% beta)) * null$score_residuals
      reference <- sum(r) / sqrt(sum((r - mean(r))^2))
      before <- null
      expect_equal(glm_score_statistic(null, exposure), reference,
                   tolerance = 1e-11)
      expect_identical(null, before)
    }
  }
})

test_that("score orientation, units, offsets and extreme scales are handled", {
  null <- .score_statistic_null()
  exposure <- seq_len(null$nobs) %% 7L
  expected <- glm_score_statistic(null, exposure)
  for (scale in c(1e-300, 1e-200, 1e200, 1e300)) {
    expect_equal(glm_score_statistic(null, scale * exposure), expected,
                 tolerance = 1e-11)
  }
  expect_equal(glm_score_statistic(null, -exposure), -expected, tolerance = 1e-11)
  expect_equal(glm_score_statistic(null, 1e15 + exposure), expected,
               tolerance = 1e-11)
  extreme <- rep(c(-1e308, 1e308), length.out = null$nobs)
  expect_equal(glm_score_statistic(null, extreme),
               glm_score_statistic(null, extreme / 1e308), tolerance = 1e-11)
})

test_that("exposure identities are checked without silently reordering", {
  null <- .score_statistic_null()
  exposure <- sin(seq_len(null$nobs))
  named <- stats::setNames(exposure, null$row_names)
  expect_identical(glm_score_statistic(null, named),
                   glm_score_statistic(null, exposure))
  expect_error(glm_score_statistic(null, rev(named)), "row order")
  names(named)[1L] <- "unknown"
  expect_error(glm_score_statistic(null, named), "row order")
})

test_that("nuisance aliases and undefined score variances fail clearly", {
  null <- .score_statistic_null()
  z <- null$projection$z[, "z"]
  for (exposure in list(rep(0, null$nobs), rep(1e200, null$nobs),
                        z, 3 + 2 * z, z * 1e-200,
                        z + 1e-10 * sin(seq_len(null$nobs)))) {
    expect_error(glm_score_statistic(null, exposure), "identifiable variation")
  }
  expect_true(is.finite(glm_score_statistic(
    null, z + 1e-4 * sin(seq_len(null$nobs))
  )))
  null$score_residuals[] <- 0
  expect_error(glm_score_statistic(null, sin(seq_len(null$nobs))),
               "positive empirical variance")
})

test_that("malformed nulls, exposures and tolerances are rejected", {
  null <- .score_statistic_null()
  n <- null$nobs
  exposure <- sin(seq_len(n))
  for (invalid in list(NULL, null$fit, unclass(null))) {
    expect_error(glm_score_statistic(invalid, exposure), "null must")
  }
  for (field in c("nobs", "row_names", "score_residuals", "projection")) {
    broken <- null
    broken[[field]] <- NULL
    expect_error(glm_score_statistic(broken, exposure), "null must")
  }
  for (invalid in list(exposure[-1L], matrix(exposure, n), rep(TRUE, n),
                       factor(rep(0:1, length.out = n)), rep(1i, n),
                       c(NA_real_, exposure[-1L]), c(Inf, exposure[-1L]))) {
    expect_error(glm_score_statistic(null, invalid), "exposure must")
  }
  for (tol in list(0, 1, -1, NA_real_, Inf, TRUE, 1i, c(0.1, 0.2))) {
    expect_error(glm_score_statistic(null, exposure, tol), "tol must")
  }
  broken <- null
  broken$row_names <- rev(broken$row_names)
  expect_error(glm_score_statistic(broken, exposure), "null must")
})
