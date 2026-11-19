.gaussian_diagnostics_from_fit <- function(fit) {
  offset <- if (is.null(fit$offset)) rep(0, length(fit$y)) else fit$offset
  gaussian_fit_diagnostics(
    fit$y, fit$fitted.values, fit$x, fit$coefficients, offset,
    fit$family$mu.eta(fit$linear.predictors)
  )
}

.gaussian_diagnostic_fit <- function(y, z, offset = rep(0, length(y)),
                                    link = "identity") {
  data <- data.frame(y, z, fixed_offset = offset)
  glm(y ~ z + offset(fixed_offset), data = data, family = gaussian(link),
      x = TRUE, y = TRUE, model = TRUE,
      control = glm.control(epsilon = 1e-10, maxit = 100))
}

test_that("Gaussian diagnostics reject nonconstant exact nuisance fits", {
  z <- seq(-1, 1, length.out = 40)
  fit <- .gaussian_diagnostic_fit(2 + 3 * z, z)
  expect_true(fit$converged)
  expect_false(fit$boundary)
  expect_true(any(fit$y != fit$y[1L]))
  expect_error(.gaussian_diagnostics_from_fit(fit), "roundoff")

  group <- factor(rep(0:3, each = 10))
  data <- data.frame(group, y = 3 + 2 * as.numeric(group))
  grouped <- glm(y ~ group, data = data, family = gaussian(),
                 x = TRUE, y = TRUE, model = TRUE)
  expect_true(grouped$converged)
  expect_false(grouped$boundary)
  expect_error(.gaussian_diagnostics_from_fit(grouped), "roundoff")
})

test_that("Gaussian exact-fit detection covers supported links and fixed offsets", {
  z <- seq(-1, 1, length.out = 40)
  fixed_offset <- 0.2 * sin(seq_along(z))
  for (link in c("identity", "log", "inverse")) {
    family <- gaussian(link)
    for (offset in list(rep(0, length(z)), fixed_offset)) {
      eta <- 2 + 0.3 * z + offset
      fit <- .gaussian_diagnostic_fit(family$linkinv(eta), z, offset, link)
      expect_true(fit$converged)
      expect_error(.gaussian_diagnostics_from_fit(fit), "roundoff")
    }
  }
})

test_that("Gaussian diagnostic decisions do not impose outcome-unit floors", {
  z <- seq(-1, 1, length.out = 80)
  noise <- 0.2 * sin(seq_along(z) * sqrt(2))
  y <- 2 + 0.7 * z + noise
  baseline <- .gaussian_diagnostics_from_fit(.gaussian_diagnostic_fit(y, z))
  for (scale in c(1e-150, 1e-50, 1, 1e50, 1e150)) {
    noisy <- .gaussian_diagnostic_fit(scale * y, z)
    diagnostic <- .gaussian_diagnostics_from_fit(noisy)
    expect_gt(diagnostic$residual_norm, diagnostic$roundoff_bound)
    expect_equal(diagnostic$residual_norm / scale, baseline$residual_norm,
                 tolerance = 1e-12)
    expect_equal(diagnostic$relative_residual_norm, baseline$relative_residual_norm,
                 tolerance = 1e-12)
    exact <- .gaussian_diagnostic_fit(scale * (2 + 0.7 * z), z)
    expect_error(.gaussian_diagnostics_from_fit(exact), "roundoff")
  }
})

test_that("resolved small noise and large common offsets remain supported", {
  z <- seq(-1, 1, length.out = 80)
  pattern <- sin(seq_along(z) * sqrt(2))
  small_noise <- .gaussian_diagnostic_fit(2 + 0.7 * z + 1e-10 * pattern, z)
  small_diagnostic <- .gaussian_diagnostics_from_fit(small_noise)
  expect_gt(small_diagnostic$residual_norm, small_diagnostic$roundoff_bound)
  for (scale in c(1e-150, 1e150)) {
    scaled <- .gaussian_diagnostic_fit(scale * (2 + 0.7 * z + 1e-10 * pattern), z)
    diagnostic <- .gaussian_diagnostics_from_fit(scaled)
    expect_gt(diagnostic$residual_norm, diagnostic$roundoff_bound)
  }

  shifted <- .gaussian_diagnostic_fit(1e12 + 2 + 0.7 * z + 0.1 * pattern, z)
  shifted_diagnostic <- .gaussian_diagnostics_from_fit(shifted)
  expect_gt(shifted_diagnostic$residual_norm, shifted_diagnostic$roundoff_bound)

  fixed_offset <- 1e12 + 0.2 * cos(seq_along(z))
  with_offset <- .gaussian_diagnostic_fit(fixed_offset + 2 + 0.7 * z + 0.1 * pattern,
                                          z, fixed_offset)
  offset_diagnostic <- .gaussian_diagnostics_from_fit(with_offset)
  expect_gt(offset_diagnostic$residual_norm, offset_diagnostic$roundoff_bound)
})

test_that("Gaussian diagnostic scale accounts for design-column rescaling", {
  z <- seq(-1, 1, length.out = 80)
  y <- 2 + 0.7 * z + 0.2 * sin(seq_along(z) * sqrt(2))
  baseline <- .gaussian_diagnostics_from_fit(.gaussian_diagnostic_fit(y, z))
  for (scale in c(1e-4, 1e4)) {
    diagnostic <- .gaussian_diagnostics_from_fit(.gaussian_diagnostic_fit(y, scale * z))
    expect_equal(diagnostic$residual_norm, baseline$residual_norm, tolerance = 1e-12)
    expect_equal(diagnostic$roundoff_bound, baseline$roundoff_bound, tolerance = 1e-12)
  }
})

test_that("Gaussian log and inverse links retain resolved residual noise", {
  z <- seq(-1, 1, length.out = 80)
  pattern <- sin(seq_along(z) * sqrt(2))
  for (link in c("log", "inverse")) {
    family <- gaussian(link)
    mean <- family$linkinv(2 + 0.3 * z)
    small_noise_fit <- .gaussian_diagnostic_fit(mean + 1e-8 * pattern, z,
                                                link = link)
    small_diagnostic <- .gaussian_diagnostics_from_fit(small_noise_fit)
    expect_gt(small_diagnostic$residual_norm, small_diagnostic$roundoff_bound)
    y <- mean + 0.02 * pattern
    for (scale in c(1e-30, 1, 1e30)) {
      fit <- .gaussian_diagnostic_fit(scale * y, z, link = link)
      expect_true(fit$converged)
      diagnostic <- .gaussian_diagnostics_from_fit(fit)
      expect_gt(diagnostic$residual_norm, diagnostic$roundoff_bound)
    }
  }
})

test_that("Gaussian score preparation enforces and records the residual check", {
  z <- seq(-1, 1, length.out = 40)
  exact <- .gaussian_diagnostic_fit(2 + 3 * z, z)
  expect_error(prepare_glm_score(exact), "roundoff")
  noisy <- .gaussian_diagnostic_fit(2 + 3 * z + 0.2 * sin(seq_along(z)), z)
  prepared <- prepare_glm_score(noisy)
  diagnostic <- prepared$diagnostics$gaussian_residual_check
  expect_named(diagnostic,
                 c("residual_norm", "roundoff_bound", "relative_residual_norm"))
  expect_gt(diagnostic$residual_norm, diagnostic$roundoff_bound)
  expect_true(all(is.finite(unlist(diagnostic))))
})
