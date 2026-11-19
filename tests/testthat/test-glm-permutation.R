.permutation_null <- function(family = gaussian()) {
  i <- seq_len(45)
  data <- data.frame(y = sin(i * 0.67) + cos(i * 0.13) + 0.2 * i / 45,
                     z = cos(i * 0.31), row.names = paste0("sample", i))
  if (family$family == "poisson") {
    data$y <- (i * 7L %% 11L) + as.integer(data$z > 0)
  }
  fit_glm_null(y ~ z, data, family = family)
}

test_that("fixed-budget calibration matches independently permuted R scores", {
  null <- .permutation_null()
  exposure <- as.numeric(seq_len(null$nobs) %% 3L == 0L)
  B <- 29L
  set.seed(309)
  r_statistics <- replicate(B, {
    x <- exposure[sample.int(null$nobs, null$nobs, replace = FALSE)]
    r <- residualize_exposure(null$projection, x) * null$score_residuals
    sum(r) / sqrt(sum((r - mean(r))^2))
  })
  reference_seed <- .Random.seed
  r <- residualize_exposure(null$projection, exposure) * null$score_residuals
  observed <- sum(r) / sqrt(sum((r - mean(r))^2))
  before <- null
  for (alternative in c("two.sided", "greater", "less")) {
    set.seed(309)
    result <- calibrate_glm_permutation(null, exposure, B = B,
      alternative = alternative, batch_size = 7L, keep_statistics = TRUE)
    expected_count <- switch(alternative,
      two.sided = sum(abs(r_statistics) >= abs(observed)),
      greater = sum(r_statistics >= observed),
      less = sum(r_statistics <= observed)
    )
    expect_equal(result$statistic, observed, tolerance = 1e-11)
    expect_equal(result$permuted_statistics, r_statistics, tolerance = 1e-11)
    expect_equal(result$p.value, (1 + expected_count) / (B + 1))
    expect_equal(result$exceedances, expected_count)
    expect_identical(result$alternative, alternative)
    expect_identical(result$B, B)
    expect_identical(.Random.seed, reference_seed)
    expect_identical(null, before)
  }
})

test_that("batch size and optional retention leave calibration and RNG unchanged", {
  null <- .permutation_null()
  exposure <- stats::setNames(sin(seq_len(null$nobs)), null$row_names)
  run <- function(batch_size, keep_statistics = TRUE) {
    set.seed(991)
    result <- calibrate_glm_permutation(null, exposure, B = 31L,
      batch_size = batch_size, keep_statistics = keep_statistics)
    list(result = result, seed = .Random.seed)
  }
  reference <- run(1L)
  expect_identical(run(7L), reference)
  expect_identical(run(100L), reference)
  streamed <- run(8L, FALSE)
  expect_null(streamed$result$permuted_statistics)
  expect_true("permuted_statistics" %in% names(streamed$result))
  reference$result$permuted_statistics <- NULL
  # Explicit NULL is retained as part of the result contract.
  streamed$result$permuted_statistics <- NULL
  expect_identical(streamed, reference)
})

test_that("identity draws, repeated permutations and ties are counted", {
  data <- data.frame(y = c(-3, 0, 2, 4), row.names = letters[1:4])
  null <- fit_glm_null(y ~ 1, data)
  exposure <- c(0, 0, 1, 1)
  set.seed(712)
  indices <- replicate(80L, sample.int(4L, 4L, replace = FALSE))
  draws <- apply(indices, 2L, function(index) exposure[index])
  identity_draws <- colSums(draws != exposure) == 0L
  expect_true(any(identity_draws))
  expect_lte(nrow(unique(t(draws))), 6L)
  set.seed(712)
  result <- calibrate_glm_permutation(null, exposure, B = 80L,
    batch_size = 9L, keep_statistics = TRUE, alternative = "greater")
  expect_true(all(result$permuted_statistics[identity_draws] == result$statistic))
  expect_equal(result$exceedances,
    sum(result$permuted_statistics >= result$statistic))
  expect_equal(result$p.value,
    permutation_p_value(result$statistic, result$permuted_statistics, "greater"))
})

test_that("nonconstant GLM weights are reused and each draw is studentized", {
  null <- .permutation_null(poisson())
  exposure <- sin(seq_len(null$nobs))
  expect_gt(diff(range(null$working_weights)), 0.1)
  set.seed(41)
  reference <- replicate(13L, {
    x <- exposure[sample.int(null$nobs, null$nobs, replace = FALSE)]
    r <- residualize_exposure(null$projection, x) * null$score_residuals
    sum(r) / sqrt(sum((r - mean(r))^2))
  })
  set.seed(41)
  result <- calibrate_glm_permutation(null, exposure, B = 13L,
    batch_size = 4L, keep_statistics = TRUE)
  expect_equal(result$permuted_statistics, reference, tolerance = 1e-11)
  one <- calibrate_glm_permutation(null, exposure, B = 1L,
    keep_statistics = TRUE)
  expect_length(one$permuted_statistics, 1L)
  expect_true(one$p.value %in% c(0.5, 1))
})

test_that("calibration uses cached null quantities without refitting", {
  null <- .permutation_null()
  exposure <- sin(seq_len(null$nobs))
  set.seed(14)
  reference <- calibrate_glm_permutation(null, exposure, B = 3L)
  # Stored fit/call are provenance; scoring needs only the prepared cache.
  null$fit <- NULL
  null$call <- quote(stop("a refit must not occur"))
  set.seed(14)
  expect_identical(calibrate_glm_permutation(null, exposure, B = 3L), reference)
})

test_that("undefined permutation scores stop calibration without redrawing", {
  data <- data.frame(y = c(-2, 1, 0, 3, 1, 4), group = rep(0:1, each = 3L))
  null <- fit_glm_null(y ~ group, data)
  exposure <- c(0, 1, 0, 1, 0, 1)
  expect_true(is.finite(glm_score_statistic(null, exposure)))
  # This seed is selected independently from labels: record the first draw
  # whose binary exposure is exactly in the nuisance span.
  set.seed(157)
  first_alias <- 0L
  for (b in seq_len(200L)) {
    x <- exposure[sample.int(6L, 6L, replace = FALSE)]
    if (all(x == data$group) || all(x == 1 - data$group)) {
      first_alias <- b
      break
    }
  }
  expect_gt(first_alias, 0L)
  expected_seed <- .Random.seed
  set.seed(157)
  expect_error(calibrate_glm_permutation(null, exposure, B = 200L, batch_size = 1L),
    paste0("Permutations ", first_alias, "-", first_alias, ":.*identifiable variation"))
  expect_identical(.Random.seed, expected_seed)
})

test_that("invalid calibration inputs fail before consuming randomness", {
  null <- .permutation_null()
  exposure <- sin(seq_len(null$nobs))
  set.seed(42)
  before <- .Random.seed
  for (bad in list(0, -1, 1.5, NA_real_, Inf, TRUE, 1i, numeric(), c(1, 2),
                   .Machine$integer.max + 1)) {
    expect_error(calibrate_glm_permutation(null, exposure, B = bad), "B must")
    expect_error(calibrate_glm_permutation(null, exposure, batch_size = bad),
                 "batch_size must")
  }
  for (bad in list(NA, 1, "yes", c(TRUE, FALSE), matrix(TRUE))) {
    expect_error(calibrate_glm_permutation(null, exposure, keep_statistics = bad),
                 "keep_statistics must")
  }
  expect_error(calibrate_glm_permutation(null, exposure, alternative = "unknown"))
  expect_error(calibrate_glm_permutation(null, rep(1, null$nobs)),
               "identifiable variation")
  named <- stats::setNames(exposure, rev(null$row_names))
  expect_error(calibrate_glm_permutation(null, named), "row order")
  expect_identical(.Random.seed, before)
})
