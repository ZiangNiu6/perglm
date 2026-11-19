.public_glm_data <- function(poisson_response = FALSE) {
  i <- seq_len(45)
  data <- data.frame(
    z = cos(i * 0.31),
    y = sin(i * 0.67) + cos(i * 0.13) + 0.2 * i / 45,
    x = as.integer(i %% 3L == 0L),
    row.names = paste0("subject", i)
  )
  if (poisson_response) {
    data$y <- (i * 7L %% 11L) + as.integer(data$z > 0)
  }
  data
}

.public_glm_reference <- function(data, family, B, seed) {
  fit <- stats::glm(y ~ z, data = data, family = family,
                    control = stats::glm.control(epsilon = 1e-10), x = TRUE)
  z <- fit$x
  mu <- fit$fitted.values
  derivative <- family$mu.eta(fit$linear.predictors)
  variance <- family$variance(mu)
  epsilon <- (data$y - mu) * derivative / variance
  weights <- derivative^2 / variance
  weighted_qr <- qr(sqrt(weights) * z)
  score <- function(exposure) {
    beta <- qr.coef(weighted_qr, sqrt(weights) * exposure)
    contributions <- (exposure - drop(z %*% beta)) * epsilon
    sum(contributions) /
      sqrt(sum((contributions - mean(contributions))^2))
  }
  observed <- score(data$x)
  set.seed(seed)
  statistics <- replicate(B, {
    score(data$x[sample.int(nrow(data), nrow(data), replace = FALSE)])
  })
  list(observed = observed, statistics = statistics, seed = .Random.seed)
}

test_that("the exported test agrees with independent GLM permutation calculations", {
  B <- 19L
  seed <- 337L
  for (family in list(stats::gaussian(), stats::poisson())) {
    dat <- .public_glm_data(family$family == "poisson")
    reference <- .public_glm_reference(dat, family, B, seed)
    for (alternative in c("two.sided", "greater", "less")) {
      expected_count <- switch(alternative,
        two.sided = sum(abs(reference$statistics) >= abs(reference$observed)),
        greater = sum(reference$statistics >= reference$observed),
        less = sum(reference$statistics <= reference$observed)
      )
      set.seed(seed)
      result <- perglm_test(y ~ z, data = dat, exposure = dat$x,
        family = family, B = B, alternative = alternative,
        batch_size = 6L, keep_statistics = TRUE,
        control = list(epsilon = 1e-10))
      expect_equal(unname(result$statistic), reference$observed,
                   tolerance = 1e-11)
      expect_equal(result$permuted_statistics, reference$statistics,
                   tolerance = 1e-11)
      expect_equal(result$exceedances, expected_count)
      expect_equal(result$p.value, (1 + expected_count) / (B + 1))
      expect_identical(result$alternative, alternative)
      expect_identical(.Random.seed, reference$seed)
    }
  }
})

test_that("the public result follows the standard htest interface", {
  dat <- .public_glm_data()
  set.seed(402)
  result <- perglm_test(y ~ z, data = dat, exposure = dat$x, B = 7L)
  expect_identical(getExportedValue("perglm", "perglm_test"), perglm_test)
  expect_identical(class(result), "htest")
  expect_identical(names(result$statistic), "T")
  expect_identical(names(result$parameter), "B")
  expect_equal(unname(result$parameter), 7L)
  expect_equal(result$B, 7L)
  expect_equal(result$nobs, nrow(dat))
  expect_identical(result$family, "gaussian")
  expect_identical(result$link, "identity")
  expect_null(result$theta)
  expect_match(tolower(result$method), "gaussian")
  expect_match(tolower(result$method), "identity")
  expect_match(result$data.name, "y ~ z", fixed = TRUE)
  expect_match(result$data.name, "dat$x", fixed = TRUE)
  expect_equal(result$null.value, c("exposure coefficient" = 0))
  expect_identical(result$alternative, "two.sided")
  expect_type(result$null_diagnostics, "list")
  expect_true("permuted_statistics" %in% names(result))
  expect_null(result$permuted_statistics)
  expect_false(any(c("fit", "null", "projection") %in% names(result)))
  output <- paste(capture.output(print(result)), collapse = "\n")
  expect_match(output, "p-value")
  expect_match(output, "T =", fixed = TRUE)
  expect_match(output, "exposure coefficient", fixed = TRUE)
})

test_that("public batching and retention preserve results and random draws", {
  dat <- .public_glm_data(TRUE)
  run <- function(batch_size, keep_statistics = TRUE) {
    set.seed(218)
    result <- perglm_test(y ~ z, data = dat, exposure = dat$x,
      family = stats::poisson(), B = 23L, batch_size = batch_size,
      keep_statistics = keep_statistics)
    list(result = result, seed = .Random.seed)
  }
  reference <- run(1L)
  expect_identical(run(7L), reference)
  expect_identical(run(100L), reference)
  discarded <- run(6L, FALSE)
  expect_true("permuted_statistics" %in% names(discarded$result))
  expect_null(discarded$result$permuted_statistics)
  reference$result$permuted_statistics <- NULL
  discarded$result$permuted_statistics <- NULL
  expect_identical(discarded, reference)
})

test_that("public exposure vectors are evaluated and checked against model rows", {
  dat <- .public_glm_data()
  named <- stats::setNames(dat$x, row.names(dat))
  set.seed(719)
  unnamed_result <- perglm_test(y ~ z, dat, dat$x, B = 5L,
                                 keep_statistics = TRUE)
  set.seed(719)
  named_result <- perglm_test(y ~ z, dat, named, B = 5L,
                               keep_statistics = TRUE)
  expect_identical(named_result$statistic, unnamed_result$statistic)
  expect_identical(named_result$p.value, unnamed_result$p.value)
  expect_identical(named_result$permuted_statistics,
                   unnamed_result$permuted_statistics)
  expect_error(perglm_test(y ~ z, dat, rev(named), B = 1L), "row order")
  expect_error(perglm_test(y ~ z, dat, "x", B = 1L), "exposure")
})

test_that("invalid public inputs and unidentified tests fail before sampling", {
  dat <- .public_glm_data()
  incomplete <- dat
  incomplete$y[2L] <- NA_real_
  missing_exposure <- dat$x
  missing_exposure[2L] <- NA_real_
  flat <- data.frame(y = rep(0, 6), x = rep(0:1, 3))
  set.seed(151)
  before <- .Random.seed
  expect_error(perglm_test(y ~ z, incomplete, incomplete$x, B = 1L))
  expect_error(perglm_test(y ~ z, dat, missing_exposure, B = 1L), "exposure")
  expect_error(perglm_test(y ~ z - 1, dat, dat$x, B = 1L), "intercept")
  expect_error(perglm_test(y ~ z, dat, dat$x, family = stats::quasipoisson(),
                            B = 1L), "family")
  expect_error(perglm_test(y ~ z, dat, dat$z, B = 1L), "identifiable variation")
  expect_error(perglm_test(y ~ z, dat, rep(1, nrow(dat)), B = 1L),
               "identifiable variation")
  expect_error(perglm_test(y ~ 1, flat, flat$x, B = 1L),
               "positive empirical variance")
  expect_error(perglm_test(y ~ z, dat, dat$x, B = 0), "B must")
  expect_error(perglm_test(y ~ z, dat, dat$x, B = 1L, keep_statistics = NA),
               "keep_statistics")
  expect_identical(.Random.seed, before)
})

test_that("negative-binomial tests report estimated and fixed size metadata", {
  set.seed(431)
  dat <- data.frame(z1 = rnorm(150), z2 = rnorm(150))
  mu <- exp(0.3 + 0.25 * dat$z1 - 0.15 * dat$z2)
  dat$y <- rnbinom(150, mu = mu, size = 1.7)
  dat$x <- rbinom(150, 1L, 0.3)
  estimated <- perglm_test(y ~ z1 + z2, dat, dat$x,
    family = "negative.binomial", B = 3L, control = list(maxit = 100L))
  fixed <- perglm_test(y ~ z1 + z2, dat, dat$x,
    family = MASS::negative.binomial(theta = 1.7, link = "log"), B = 3L)
  for (result in list(estimated, fixed)) {
    expect_s3_class(result, "htest")
    expect_identical(result$family, "negative.binomial")
    expect_identical(result$link, "log")
    expect_length(result$theta, 1L)
    expect_true(is.finite(result$theta))
    expect_gt(result$theta, 0)
    expect_true(is.finite(result$p.value))
  }
  expect_identical(fixed$theta, 1.7)
})

test_that("a public family factory gives the same test as its family object", {
  dat <- .public_glm_data(TRUE)
  factory <- function() stats::poisson(link = "log")
  set.seed(829)
  from_factory <- perglm_test(y ~ z, dat, dat$x, family = factory, B = 3L)
  set.seed(829)
  from_object <- perglm_test(y ~ z, dat, dat$x,
    family = stats::poisson(link = "log"), B = 3L)
  expect_identical(from_factory, from_object)
})
