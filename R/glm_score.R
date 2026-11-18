#' Prepare score quantities from a fitted null GLM
#'
#' Use a null model containing the nuisance covariates only. The fit must retain
#' its model frame, design and response (`model = TRUE, x = TRUE, y = TRUE`).
#' Observations must be complete and unweighted, with an intercept and positive
#' residual degrees of freedom. Fixed offsets are supported. Binomial responses
#' must be Bernoulli; grouped counts and proportions are not supported.
#'
#' Standard Gaussian, binomial, Poisson, Gamma, inverse Gaussian and MASS
#' negative-binomial families are supported with their standard named links.
#' Common dispersion factors are omitted because they cancel from empirical
#' studentization. Negative-binomial size is retained because it changes
#' relative working weights. All quantities are computed at the final means.
#'
#' The returned diagnostics include the maximum normalized nuisance score,
#' `abs(crossprod(z[, j], epsilon)) / (norm(z[, j]) * norm(epsilon))`.
#' This is recorded rather than used as an additional convergence cutoff.
#' Preparation does not inspect an exposure. The testing layer must separately
#' check that the exposure has residual variation after nuisance adjustment.
#'
#' @param fit A standard `stats::glm()` or `MASS::glm.nb()` null-model fit.
#' @param tol Relative rank tolerance passed to `prepare_projection()`.
#' @return An internal `perglm_glm_score` object holding the fit, working family,
#'   score residuals, working weights, projection, row identities and diagnostics.
#' @noRd
prepare_glm_score <- function(fit, tol = 1e-7) {
  expected_class <- if (inherits(fit, "negbin")) {
    c("negbin", "glm", "lm")
  } else {
    c("glm", "lm")
  }
  if (!is.list(fit) || !identical(class(fit), expected_class) ||
      !identical(fit$method, "glm.fit")) {
    stop("fit must be a standard stats::glm() or MASS::glm.nb() result.",
         call. = FALSE)
  }
  if (!isTRUE(fit$converged) || isTRUE(fit$boundary)) {
    stop("The null fit must have converged to an interior solution.", call. = FALSE)
  }
  if (length(fit$th.warn)) {
    stop("The negative-binomial size estimate has a convergence warning.",
         call. = FALSE)
  }
  if (!is.data.frame(fit$model) || is.null(fit$x) || is.null(fit$y)) {
    stop("Refit the null model with model = TRUE, x = TRUE, and y = TRUE.",
         call. = FALSE)
  }
  if (!is.null(fit$na.action) || !is.null(attr(fit$model, "na.action")) ||
      !is.null(fit$call$subset)) {
    stop("The null fit must use complete observations without omitted or subset rows.",
         call. = FALSE)
  }

  z <- fit$x
  if (!is.matrix(z) || !is.numeric(z) || is.complex(z) ||
      any(!is.finite(z)) || nrow(z) < 2L || ncol(z) < 1L) {
    stop("The retained design must be a finite real numeric matrix.", call. = FALSE)
  }
  n <- nrow(z)
  p <- ncol(z)
  if (!identical(attr(fit$terms, "intercept"), 1L) ||
      !"(Intercept)" %in% colnames(z) ||
      any(z[, "(Intercept)"] != 1)) {
    stop("The null design must include an intercept.", call. = FALSE)
  }
  if (length(fit$rank) != 1L || !is.finite(fit$rank) ||
      fit$rank != p || n <= p) {
    stop("The null design must have full column rank and positive residual degrees of freedom.",
         call. = FALSE)
  }
  if (!glm_numeric_vector(fit$prior.weights, n) ||
      any(fit$prior.weights != 1)) {
    stop("Only unit prior weights and Bernoulli trial counts are supported.",
         call. = FALSE)
  }
  rows <- row.names(fit$model)
  if (nrow(fit$model) != n || !identical(rownames(z), rows)) {
    stop("The retained null-model rows are not aligned.", call. = FALSE)
  }
  for (value in list(fit$y, fit$fitted.values, fit$linear.predictors)) {
    if (!glm_numeric_vector(value, n)) {
      stop("The retained response, fitted means, and predictors must be finite real vectors.",
           call. = FALSE)
    }
    if (!is.null(names(value)) && !identical(names(value), rows)) {
      stop("The retained null-model rows are not aligned.", call. = FALSE)
    }
  }
  if (!glm_numeric_vector(fit$coefficients, p) ||
      (!is.null(names(fit$coefficients)) &&
       !identical(names(fit$coefficients), colnames(z)))) {
    stop("The null coefficients must be finite and aligned with the design.",
         call. = FALSE)
  }
  offset <- if (is.null(fit$offset)) rep(0, n) else fit$offset
  if (!glm_numeric_vector(offset, n)) {
    stop("The fixed offset must be a finite real vector aligned with the design.",
         call. = FALSE)
  }

  if (inherits(fit, "negbin") &&
      (!glm_numeric_vector(fit$theta, 1L) || fit$theta <= 0)) {
    stop("The negative-binomial size theta must be finite and strictly positive.",
         call. = FALSE)
  }
  working <- standard_glm_family(fit$family,
    theta = if (inherits(fit, "negbin")) fit$theta else NULL
  )
  if (inherits(fit, "negbin") && working$name != "negative.binomial") {
    stop("A negbin fit must have a negative-binomial family.", call. = FALSE)
  }
  family <- working$family
  y <- fit$y
  mu <- fit$fitted.values
  eta <- fit$linear.predictors
  if (working$name == "binomial" && any(!y %in% c(0, 1))) {
    stop("Binomial responses must be Bernoulli observations coded as zero or one.",
         call. = FALSE)
  }
  if (working$name %in% c("poisson", "negative.binomial") &&
      any(y < 0 | y != floor(y))) {
    stop("Poisson and negative-binomial responses must be nonnegative integer counts.",
         call. = FALSE)
  }
  if (working$name %in% c("Gamma", "inverse.gaussian") && any(y <= 0)) {
    stop("Gamma and inverse-Gaussian responses must be strictly positive.",
         call. = FALSE)
  }
  if (!isTRUE(family$validmu(mu)) || !isTRUE(family$valideta(eta))) {
    stop("The fitted means and predictors must be inside the working family domain.",
         call. = FALSE)
  }
  expected_mu <- family$linkinv(eta)
  expected_eta <- drop(z %*% fit$coefficients) + offset
  if (!glm_vectors_close(mu, expected_mu, 1e-10) ||
      !glm_vectors_close(eta, expected_eta, 1e-8)) {
    stop("The fitted means, predictors, coefficients, and offset are inconsistent.",
         call. = FALSE)
  }

  derivative <- family$mu.eta(eta)
  variance <- family$variance(mu)
  if (!glm_numeric_vector(derivative, n) || any(derivative == 0) ||
      !glm_numeric_vector(variance, n) || any(variance <= 0)) {
    stop("Working link derivatives must be finite and nonzero, and variances positive.",
         call. = FALSE)
  }
  score_residuals <- (y - mu) * (derivative / variance)
  working_weights <- (derivative / sqrt(variance))^2
  if (any(!is.finite(score_residuals)) || any(!is.finite(working_weights)) ||
      any(working_weights <= 0)) {
    stop("The working score residuals and weights must be finite, with positive weights.",
         call. = FALSE)
  }
  projection <- prepare_projection(z, working_weights, tol = tol)

  # Normalize both vectors before inner products to avoid overflow in this
  # diagnostic; it does not replace the actual signed score contributions.
  epsilon_scale <- max(abs(score_residuals))
  score_error <- 0
  if (epsilon_scale > 0) {
    scaled_epsilon <- score_residuals / epsilon_scale
    scaled_z <- sweep(z, 2L, apply(abs(z), 2L, max), "/")
    score_error <- max(abs(drop(crossprod(scaled_z, scaled_epsilon))) /
      (sqrt(colSums(scaled_z^2)) * sqrt(sum(scaled_epsilon^2))))
  }

  structure(list(
    fit = fit,
    projection = projection,
    score_residuals = score_residuals,
    working_weights = working_weights,
    family = family,
    family_name = working$name,
    theta = working$theta,
    nobs = n,
    row_names = rows,
    diagnostics = list(converged = TRUE, iterations = fit$iter,
                       max_score_correlation = score_error)
  ), class = "perglm_glm_score")
}

#' @noRd
glm_numeric_vector <- function(value, n) {
  is.numeric(value) && !is.complex(value) && is.null(dim(value)) &&
    length(value) == n && all(is.finite(value))
}

#' @noRd
glm_vectors_close <- function(x, y, tolerance) {
  glm_numeric_vector(y, length(x)) &&
    all(abs(x - y) <= tolerance * pmax(1, abs(x), abs(y)))
}

#' Reconstruct a supported working family
#'
#' Rebuilding the family avoids stale negative-binomial shape closures. Shape
#' is read directly from MASS's family closure rather than its rounded label.
#' @noRd
standard_glm_family <- function(family, theta = NULL) {
  if (!inherits(family, "family") || !is.character(family$family) ||
      length(family$family) != 1L || is.na(family$family) ||
      !is.character(family$link) || length(family$link) != 1L || is.na(family$link)) {
    stop("A supported standard GLM family with a named link is required.",
         call. = FALSE)
  }
  name <- family$family
  links <- list(
    gaussian = c("identity", "log", "inverse"),
    binomial = c("logit", "probit", "cauchit", "log", "identity", "cloglog"),
    poisson = c("log", "identity", "sqrt"),
    Gamma = c("inverse", "identity", "log"),
    inverse.gaussian = c("1/mu^2", "inverse", "identity", "log"),
    negative.binomial = c("log", "identity", "sqrt")
  )
  if (startsWith(name, "Negative Binomial(")) {
    name <- "negative.binomial"
    if (is.null(theta)) {
      variance_environment <- if (is.function(family$variance)) {
        environment(family$variance)
      } else NULL
      if (!is.environment(variance_environment) ||
          !exists(".Theta", envir = variance_environment, inherits = FALSE)) {
        stop("Use MASS::negative.binomial() to specify a fixed size.", call. = FALSE)
      }
      theta <- get(".Theta", envir = variance_environment, inherits = FALSE)
    }
    if (!glm_numeric_vector(theta, 1L) || theta <= 0) {
      stop("The negative-binomial size theta must be finite and strictly positive.",
           call. = FALSE)
    }
  }
  if (!name %in% names(links) || !family$link %in% links[[name]]) {
    stop("Unsupported GLM family or link; use a standard Gaussian, binomial, Poisson, Gamma, inverse-Gaussian, or negative-binomial family.",
         call. = FALSE)
  }
  constructor <- switch(name,
    gaussian = stats::gaussian,
    binomial = stats::binomial,
    poisson = stats::poisson,
    Gamma = stats::Gamma,
    inverse.gaussian = stats::inverse.gaussian,
    negative.binomial = MASS::negative.binomial
  )
  rebuilt <- if (name == "negative.binomial") {
    constructor(theta = theta, link = family$link)
  } else constructor(link = family$link)
  list(family = rebuilt, name = name, theta = theta)
}
