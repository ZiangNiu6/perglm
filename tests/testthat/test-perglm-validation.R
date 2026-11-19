.validation_positive_data <- function() {
  i <- seq_len(60)
  group <- factor(rep(c("a", "b", "c"), 20))
  z <- seq(-0.8, 0.8, length.out = length(i))
  data.frame(
    y = 2 + 0.2 * z + c(-0.1, 0, 0.1)[as.integer(group)] +
      0.25 * sin(i * 1.13) + 0.15 * cos(i * 0.71),
    z = z, group = group, off = 0.03 * sin(i * 0.43),
    x = as.integer(i %% 4L == 0L),
    row.names = paste0("unit", i)
  )
}

.validation_family_cases <- function() {
  positive <- .validation_positive_data()
  bernoulli <- positive
  bernoulli$y <- factor(ifelse(seq_len(nrow(bernoulli)) %% 5L %in% c(0L, 2L),
                               "yes", "no"), levels = c("no", "yes"))
  set.seed(431)
  negative_binomial <- data.frame(
    z1 = rnorm(150), z2 = rnorm(150),
    group = factor(rep(c("a", "b", "c"), 50)),
    off = 0.1 * sin(seq_len(150) * 0.61),
    row.names = paste0("nb", seq_len(150))
  )
  mu <- exp(0.3 + 0.25 * negative_binomial$z1 -
              0.15 * negative_binomial$z2 + negative_binomial$off)
  negative_binomial$y <- rnbinom(150, mu = mu, size = 1.7)
  negative_binomial$x <- rbinom(150, 1L, 0.3)
  positive_formula <- y ~ z + group + offset(off)
  nb_formula <- y ~ z1 + z2 + group + offset(off)
  list(
    binomial_factor_logit = list(data = bernoulli, formula = positive_formula,
      family = stats::binomial(), alternative = "two.sided"),
    Gamma_inverse = list(data = positive, formula = positive_formula,
      family = stats::Gamma(), alternative = "greater", decreasing = TRUE),
    Gamma_log = list(data = positive, formula = positive_formula,
      family = stats::Gamma(link = "log"), alternative = "less"),
    inverse_Gaussian_canonical = list(data = positive, formula = positive_formula,
      family = stats::inverse.gaussian(), alternative = "greater", decreasing = TRUE),
    inverse_Gaussian_log = list(data = positive, formula = positive_formula,
      family = stats::inverse.gaussian(link = "log"), alternative = "two.sided"),
    Gaussian_inverse = list(data = positive, formula = positive_formula,
      family = stats::gaussian(link = "inverse"), alternative = "greater",
      decreasing = TRUE),
    negative_binomial_estimated = list(data = negative_binomial, formula = nb_formula,
      family = "negative.binomial", alternative = "less"),
    negative_binomial_fixed = list(data = negative_binomial, formula = nb_formula,
      family = MASS::negative.binomial(theta = 1.7), alternative = "two.sided")
  )
}

.validation_glm_reference <- function(case, B, seed) {
  control <- stats::glm.control(epsilon = 1e-10, maxit = 100L)
  estimated_nb <- identical(case$family, "negative.binomial")
  fit <- if (estimated_nb) {
    MASS::glm.nb(case$formula, data = case$data, control = control,
                  x = TRUE, y = TRUE, model = TRUE, na.action = stats::na.fail)
  } else {
    stats::glm(case$formula, data = case$data, family = case$family,
      control = control, x = TRUE, y = TRUE, model = TRUE,
      na.action = stats::na.fail, singular.ok = FALSE)
  }
  mu <- fit$fitted.values
  if (estimated_nb) {
    # Use the final size directly, not glm.nb's last fitted family closure.
    derivative <- mu
    variance <- mu + mu^2 / fit$theta
  } else {
    derivative <- case$family$mu.eta(fit$linear.predictors)
    variance <- case$family$variance(mu)
  }
  epsilon <- (fit$y - mu) * derivative / variance
  weights <- derivative^2 / variance
  score <- function(exposure, residuals = epsilon) {
    # An independent weighted least-squares fit for each exposure also checks
    # that the package's cached projection is valid for every relabeling.
    adjusted <- stats::lm.wfit(fit$x, exposure, weights)$residuals
    contributions <- adjusted * residuals
    sum_squares <- sum((contributions - mean(contributions))^2)
    stopifnot(sum_squares > 1e-18)
    sum(contributions) / sqrt(sum_squares)
  }
  observed <- score(case$data$x)
  opposite_derivative_score <- if (all(derivative < 0)) {
    score(case$data$x, residuals = -epsilon)
  } else NULL
  set.seed(seed)
  statistics <- replicate(B, {
    score(case$data$x[sample.int(nrow(case$data), nrow(case$data), replace = FALSE)])
  })
  list(observed = observed, statistics = statistics, seed = .Random.seed,
       theta = if (estimated_nb) fit$theta else NULL,
       opposite_derivative_score = opposite_derivative_score)
}

test_that("remaining families agree with independently fitted permutation scores", {
  cases <- .validation_family_cases()
  B <- 13L
  seed <- 547L
  for (name in names(cases)) {
    case <- cases[[name]]
    reference <- .validation_glm_reference(case, B, seed)
    set.seed(seed)
    result <- perglm_test(case$formula, case$data, case$data$x,
      family = case$family, B = B, alternative = case$alternative,
      batch_size = 5L, keep_statistics = TRUE,
      control = list(epsilon = 1e-10, maxit = 100L))
    expected_count <- switch(case$alternative,
      two.sided = sum(abs(reference$statistics) >= abs(reference$observed)),
      greater = sum(reference$statistics >= reference$observed),
      less = sum(reference$statistics <= reference$observed)
    )
    expect_equal(unname(result$statistic), reference$observed,
                 tolerance = 1e-10, info = name)
    expect_equal(result$permuted_statistics, reference$statistics,
                 tolerance = 1e-10, info = name)
    expect_equal(result$exceedances, expected_count, info = name)
    expect_equal(result$p.value, (1 + expected_count) / (B + 1), info = name)
    expect_identical(.Random.seed, reference$seed, info = name)
    expect_equal(result$nobs, nrow(case$data), info = name)
    if (isTRUE(case$decreasing)) {
      expect_gt(abs(reference$observed), 1e-4, label = name)
      expect_equal(unname(result$statistic), -reference$opposite_derivative_score,
                   tolerance = 1e-10, info = name)
    }
    if (identical(case$family, "negative.binomial")) {
      expect_equal(result$theta, reference$theta, tolerance = 1e-12, info = name)
    }
  }
})

test_that("exact boundary responses cannot produce a public inference result", {
  dat <- data.frame(y = rep(0, 20), z = seq(-1, 1, length.out = 20),
                     x = sin(seq_len(20)))
  set.seed(17)
  before <- .Random.seed
  for (family in list(stats::poisson(), stats::binomial(),
                      MASS::negative.binomial(theta = 1.7))) {
    expect_error(perglm_test(y ~ z, dat, dat$x, family = family, B = 9L))
  }
  dat$y[] <- 1
  expect_error(perglm_test(y ~ z, dat, dat$x, family = stats::binomial(), B = 9L))
  expect_identical(.Random.seed, before)
})

test_that("constant outcomes with a constant offset have no empirical score variance", {
  dat <- data.frame(y = rep(3, 20), z = seq(-1, 1, length.out = 20),
                     x = sin(seq_len(20)), off = rep(0.25, 20))
  set.seed(239)
  before <- .Random.seed
  expect_error(perglm_test(y ~ z, dat, dat$x, B = 3L),
               "positive empirical variance")
  expect_error(perglm_test(y ~ z + offset(off), dat, dat$x, B = 3L),
               "positive empirical variance")
  dat$y[] <- 1
  expect_error(perglm_test(y ~ z, dat, dat$x, family = stats::poisson(), B = 3L),
               "positive empirical variance")
  dat$y[] <- 2
  expect_error(perglm_test(y ~ z, dat, dat$x,
    family = MASS::negative.binomial(theta = 1.7), B = 3L),
    "positive empirical variance")
  expect_identical(.Random.seed, before)
})

test_that("a varying fixed offset can make constant outcomes informative", {
  i <- seq_len(20)
  dat <- data.frame(y = rep(3, 20), z = seq(-1, 1, length.out = 20),
                     x = cos(i * 0.43), off = 0.25 * sin(i * 0.71))
  case <- list(data = dat, formula = y ~ z + offset(off),
               family = stats::gaussian())
  reference <- .validation_glm_reference(case, B = 9L, seed = 17L)
  set.seed(17)
  result <- perglm_test(y ~ z + offset(off), dat, dat$x, B = 9L,
    batch_size = 3L, keep_statistics = TRUE,
    control = list(epsilon = 1e-10, maxit = 100L))
  expect_true(is.finite(result$statistic))
  expect_equal(unname(result$statistic), reference$observed, tolerance = 1e-11)
  expect_equal(result$permuted_statistics, reference$statistics, tolerance = 1e-11)
  expect_equal(result$p.value,
    (1 + sum(abs(reference$statistics) >= abs(reference$observed))) / 10)
  expect_identical(.Random.seed, reference$seed)
})
