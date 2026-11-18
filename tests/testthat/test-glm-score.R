.glm_score_fixture <- function(family_name = "gaussian", link = NULL,
                               seed = 912L, n = 120L) {
  set.seed(seed)
  z1 <- seq(-1, 1, length.out = n)
  z2 <- cos(seq_len(n) * 0.43)
  linear <- 0.3 + 0.25 * z1 - 0.15 * z2
  mu <- exp(linear)
  if (family_name == "Gamma" && identical(link, "inverse")) {
    mu <- 1 / (0.8 + 0.1 * z1 - 0.07 * z2)
  }
  if (family_name == "inverse.gaussian" && is.null(link)) {
    mu <- 1 / sqrt(0.7 + 0.1 * z1 - 0.05 * z2)
  }
  y <- switch(family_name,
    gaussian = 0.3 + 0.25 * z1 - 0.15 * z2 + rnorm(n, sd = 0.6),
    binomial = rbinom(n, 1, plogis(linear)),
    poisson = rpois(n, mu),
    Gamma = rgamma(n, shape = 5, scale = mu / 5),
    inverse.gaussian = {
      # Michael--Schucany--Haas construction, using only base/recommended R.
      shape <- 5
      v <- rnorm(n)^2
      candidate <- mu + mu^2 * v / (2 * shape) -
        mu / (2 * shape) * sqrt(4 * mu * shape * v + mu^2 * v^2)
      ifelse(runif(n) <= mu / (mu + candidate), candidate, mu^2 / candidate)
    },
    negative.binomial = rnbinom(n, mu = mu, size = 1.7)
  )
  data <- data.frame(y, z1, z2,
                     row.names = sprintf("observation_%03d", seq_len(n)))
  family <- switch(family_name,
    gaussian = if (is.null(link)) gaussian() else gaussian(link),
    binomial = if (is.null(link)) binomial() else binomial(link),
    poisson = if (is.null(link)) poisson() else poisson(link),
    Gamma = if (is.null(link)) Gamma() else Gamma(link),
    inverse.gaussian = if (is.null(link)) inverse.gaussian() else inverse.gaussian(link),
    negative.binomial = MASS::negative.binomial(1.7, link = "log")
  )
  glm(y ~ z1 + z2, data = data, family = family,
      x = TRUE, y = TRUE, model = TRUE,
      control = glm.control(epsilon = 1e-10, maxit = 100))
}

test_that("the adapter returns a reusable score object without changing its fit", {
  fit <- .glm_score_fixture("poisson")
  original <- fit
  before <- .Random.seed
  prepared <- prepare_glm_score(fit)
  expect_s3_class(prepared, "perglm_glm_score")
  expect_s3_class(prepared$projection, "perglm_projection")
  expect_identical(prepared$fit, original)
  expect_identical(fit, original)
  expect_identical(.Random.seed, before)
  expect_equal(prepared$nobs, nrow(fit$x))
  expect_identical(prepared$row_names, rownames(fit$model))
  expect_identical(prepared$family_name, "poisson")
  expect_identical(prepared$family$family, "poisson")
  expect_identical(prepared$family$link, "log")
  expect_null(prepared$theta)
  expect_true(prepared$diagnostics$converged)
  expect_equal(prepared$diagnostics$iterations, fit$iter)
  expect_true(is.finite(prepared$diagnostics$max_score_correlation))
  expect_gte(prepared$diagnostics$max_score_correlation, 0)
  expect_lte(prepared$diagnostics$max_score_correlation, 1 + 1e-12)
  exposure <- sin(seq_len(nrow(fit$x)))
  reference <- prepare_projection(fit$x, prepared$working_weights)
  expect_equal(residualize_exposure(prepared$projection, exposure),
               residualize_exposure(reference, exposure), tolerance = 1e-12)
})

test_that("standard families satisfy their analytical score and weight identities", {
  cases <- list(
    list(name = "gaussian", link = NULL,
         score = function(y, mu) y - mu, weight = function(mu) rep(1, length(mu))),
    list(name = "binomial", link = "logit",
         score = function(y, mu) y - mu, weight = function(mu) mu * (1 - mu)),
    list(name = "poisson", link = "log",
         score = function(y, mu) y - mu, weight = function(mu) mu),
    list(name = "Gamma", link = "log",
         score = function(y, mu) (y - mu) / mu,
         weight = function(mu) rep(1, length(mu))),
    list(name = "inverse.gaussian", link = "log",
         score = function(y, mu) (y - mu) / mu^2, weight = function(mu) 1 / mu)
  )
  for (case in cases) {
    fit <- .glm_score_fixture(case$name, case$link)
    prepared <- prepare_glm_score(fit)
    mu <- fit$fitted.values
    expect_identical(prepared$family_name, case$name)
    expect_equal(unname(prepared$score_residuals),
                 unname(case$score(fit$y, mu)), tolerance = 1e-12)
    expect_equal(unname(prepared$working_weights),
                 unname(case$weight(mu)), tolerance = 1e-12)
    expect_null(prepared$theta)
  }
})

test_that("decreasing standard links preserve the score sign and positive weights", {
  gamma_fit <- .glm_score_fixture("Gamma", "inverse")
  gamma_score <- prepare_glm_score(gamma_fit)
  expect_equal(unname(gamma_score$score_residuals),
               unname(-(gamma_fit$y - gamma_fit$fitted.values)), tolerance = 1e-12)
  expect_equal(unname(gamma_score$working_weights),
               unname(gamma_fit$fitted.values^2), tolerance = 1e-12)

  inverse_fit <- .glm_score_fixture("inverse.gaussian")
  inverse_score <- prepare_glm_score(inverse_fit)
  expect_identical(inverse_score$family$link, "1/mu^2")
  expect_equal(unname(inverse_score$score_residuals),
               unname(-0.5 * (inverse_fit$y - inverse_fit$fitted.values)),
               tolerance = 1e-12)
  expect_equal(unname(inverse_score$working_weights),
               unname(0.25 * inverse_fit$fitted.values^3), tolerance = 1e-12)
  expect_true(all(gamma_score$working_weights > 0))
  expect_true(all(inverse_score$working_weights > 0))
})

test_that("official noncanonical binary links are accepted", {
  for (link in c("probit", "cloglog", "cauchit")) {
    prepared <- prepare_glm_score(.glm_score_fixture("binomial", link))
    expect_identical(prepared$family$link, link)
    expect_true(all(is.finite(prepared$score_residuals)))
    expect_true(all(is.finite(prepared$working_weights)))
    expect_true(all(prepared$working_weights > 0))
  }
})

test_that("final fitted values determine working quantities", {
  fit <- .glm_score_fixture("poisson")
  expected <- prepare_glm_score(fit)
  fit$weights <- rep(NA_real_, length(fit$y))
  fit$family$variance <- function(mu) rep(999, length(mu))
  rebuilt <- prepare_glm_score(fit)
  expect_equal(rebuilt$score_residuals, expected$score_residuals)
  expect_equal(rebuilt$working_weights, expected$working_weights)
  expect_equal(rebuilt$family$variance(c(1, 2)), c(1, 2))
})

test_that("fixed and estimated negative-binomial shapes use the correct family", {
  skip_if_not_installed("MASS")
  fixed <- .glm_score_fixture("negative.binomial", n = 160L)
  estimated <- MASS::glm.nb(y ~ z1 + z2, data = fixed$model,
                            x = TRUE, y = TRUE, model = TRUE,
                            control = glm.control(epsilon = 1e-10, maxit = 100))
  expect_null(estimated$th.warn)
  for (fit in list(fixed, estimated)) {
    theta <- if (inherits(fit, "negbin")) fit$theta else 1.7
    prepared <- prepare_glm_score(fit)
    mu <- fit$fitted.values
    expect_identical(prepared$family_name, "negative.binomial")
    expect_equal(prepared$theta, theta)
    expect_equal(unname(prepared$score_residuals),
                 unname((fit$y - mu) / (1 + mu / theta)), tolerance = 1e-12)
    expect_equal(unname(prepared$working_weights),
                 unname(mu / (1 + mu / theta)), tolerance = 1e-12)
    expect_equal(prepared$family$variance(mu), mu + mu^2 / theta,
                 tolerance = 1e-12)
  }

  # A glm.nb object's final theta is authoritative over stale family closures.
  estimated$family <- MASS::negative.binomial(3 * estimated$theta, link = "log")
  rebuilt <- prepare_glm_score(estimated)
  expect_equal(rebuilt$theta, estimated$theta)
  expect_equal(rebuilt$family$variance(2), 2 + 4 / estimated$theta)

  # A fixed-theta glm obtains theta from the closure, not its printable label.
  fixed$family$family <- "Negative Binomial(999)"
  expect_equal(prepare_glm_score(fixed)$theta, 1.7)
})

test_that("negative-binomial shape problems are rejected", {
  skip_if_not_installed("MASS")
  fixed <- .glm_score_fixture("negative.binomial", n = 160L)
  estimated <- MASS::glm.nb(y ~ z1 + z2, data = fixed$model,
                            x = TRUE, y = TRUE, model = TRUE,
                            control = glm.control(epsilon = 1e-10, maxit = 100))
  for (theta in list(NULL, 0, -1, NA_real_, Inf, c(1, 2))) {
    bad <- estimated
    bad$theta <- theta
    expect_error(prepare_glm_score(bad))
  }
  bad <- estimated
  bad$th.warn <- "iteration limit reached"
  expect_error(prepare_glm_score(bad))
  bad <- fixed
  bad$family$variance <- function(mu) mu + mu^2 / 1.7
  environment(bad$family$variance) <- new.env(parent = baseenv())
  expect_error(prepare_glm_score(bad))
  for (theta in list(0, NA_real_, Inf, c(1, 2))) {
    bad <- fixed
    private <- new.env(parent = baseenv())
    private$.Theta <- theta
    bad$family$variance <- function(mu) mu + mu^2 / .Theta
    environment(bad$family$variance) <- private
    expect_error(prepare_glm_score(bad))
  }
})

test_that("complete row identities and fixed offsets are preserved", {
  fit <- .glm_score_fixture("poisson")
  data <- fit$model
  data$fixed_offset <- 0.15 * sin(seq_len(nrow(data)))
  offset_fit <- glm(y ~ z1 + z2 + offset(fixed_offset), data = data,
                     family = poisson(), x = TRUE, y = TRUE, model = TRUE,
                     control = glm.control(epsilon = 1e-10))
  prepared <- prepare_glm_score(offset_fit)
  expect_identical(prepared$row_names, rownames(data))
  expect_equal(unname(prepared$score_residuals),
               unname(offset_fit$y - offset_fit$fitted.values), tolerance = 1e-12)
  bad <- offset_fit
  bad$offset <- bad$offset + 0.1
  expect_error(prepare_glm_score(bad))
  bad <- offset_fit
  bad$offset[1L] <- Inf
  expect_error(prepare_glm_score(bad))

  for (field in c("y", "fitted.values", "linear.predictors")) {
    bad <- fit
    names(bad[[field]]) <- rev(rownames(fit$model))
    expect_error(prepare_glm_score(bad))
  }
  bad <- fit
  rownames(bad$x) <- rev(rownames(fit$model))
  expect_error(prepare_glm_score(bad))

  # Vector names are optional; the stored frame/design still identify rows.
  unnamed <- fit
  for (field in c("y", "fitted.values", "linear.predictors")) {
    names(unnamed[[field]]) <- NULL
  }
  expect_identical(prepare_glm_score(unnamed)$row_names, rownames(fit$model))
})

test_that("omitted rows, subsets and nonunit prior weights are rejected", {
  fit <- .glm_score_fixture("poisson")
  data <- fit$model
  subset_fit <- glm(y ~ z1 + z2, data = data, family = poisson(),
                     subset = seq_len(nrow(data)) <= 100,
                     x = TRUE, y = TRUE, model = TRUE)
  expect_error(prepare_glm_score(subset_fit))
  data$y[1L] <- NA_real_
  omitted <- glm(y ~ z1 + z2, data = data, family = poisson(),
                  na.action = na.omit, x = TRUE, y = TRUE, model = TRUE)
  expect_error(prepare_glm_score(omitted))
  bad <- fit
  attr(bad$model, "na.action") <- structure(1L, class = "omit")
  expect_error(prepare_glm_score(bad))
  for (weight in c(0, 0.5, 2, NA_real_, Inf)) {
    bad <- fit
    bad$prior.weights[1L] <- weight
    expect_error(prepare_glm_score(bad))
  }
  bad <- fit
  bad$prior.weights <- NULL
  expect_error(prepare_glm_score(bad))
})

test_that("required fit metadata, rank and intercept are validated", {
  fit <- .glm_score_fixture()
  for (field in c("x", "y", "model")) {
    bad <- fit
    bad[[field]] <- NULL
    expect_error(prepare_glm_score(bad))
  }
  for (mutate in list(
    function(x) { class(x) <- c("custom_glm", class(x)); x },
    function(x) { x$method <- "custom_fit"; x },
    function(x) { x$converged <- FALSE; x },
    function(x) { x$boundary <- TRUE; x },
    function(x) { x$rank <- x$rank - 1L; x }
  )) {
    expect_error(prepare_glm_score(mutate(fit)))
  }
  no_intercept <- glm(y ~ z1 + z2 - 1, data = fit$model,
                       x = TRUE, y = TRUE, model = TRUE)
  expect_error(prepare_glm_score(no_intercept))
  aliased <- glm(y ~ z1 + z2 + I(2 * z1), data = fit$model,
                  x = TRUE, y = TRUE, model = TRUE)
  expect_error(prepare_glm_score(aliased))
  saturated <- glm(y ~ z1 + z2, data = fit$model[1:3, ],
                    x = TRUE, y = TRUE, model = TRUE)
  expect_error(prepare_glm_score(saturated))
  for (tol in list(0, 1, NA_real_, Inf, c(0.1, 0.2), 1i)) {
    expect_error(prepare_glm_score(fit, tol = tol))
  }
})

test_that("inconsistent or nonfinite fitted quantities are rejected", {
  fit <- .glm_score_fixture("poisson")
  for (field in c("y", "fitted.values", "linear.predictors", "coefficients")) {
    for (value in c(NA_real_, Inf)) {
      bad <- fit
      bad[[field]][1L] <- value
      expect_error(prepare_glm_score(bad))
    }
  }
  bad <- fit
  bad$fitted.values[1L] <- 1.1 * bad$fitted.values[1L]
  expect_error(prepare_glm_score(bad))
  bad <- fit
  bad$linear.predictors[1L] <- bad$linear.predictors[1L] + 0.1
  expect_error(prepare_glm_score(bad))
  bad <- fit
  bad$coefficients[1L] <- bad$coefficients[1L] + 0.1
  expect_error(prepare_glm_score(bad))
})

test_that("unsupported families, custom links and invalid outcomes are rejected", {
  fit <- .glm_score_fixture("poisson")
  bad <- fit
  bad$family <- quasipoisson()
  expect_error(prepare_glm_score(bad))
  bad <- fit
  bad$family$family <- "custom_family"
  expect_error(prepare_glm_score(bad))
  bad <- fit
  bad$family$link <- "custom_log"
  expect_error(prepare_glm_score(bad))
  for (case in list(
    list(name = "binomial", value = 0.25),
    list(name = "binomial", value = 2),
    list(name = "poisson", value = -1),
    list(name = "poisson", value = 0.5),
    list(name = "Gamma", value = 0),
    list(name = "inverse.gaussian", value = -1)
  )) {
    bad <- .glm_score_fixture(case$name)
    bad$y[1L] <- case$value
    expect_error(prepare_glm_score(bad))
  }
})
