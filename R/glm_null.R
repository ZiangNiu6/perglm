#' Fit and prepare a null generalized linear model
#'
#' The formula must contain the outcome and nuisance covariates only; exclude
#' the exposure being tested. All model observations must be complete. This
#' wrapper retains the model frame, design and response and never drops rows.
#' Offsets can be included with `offset()` in the formula.
#'
#' Standard families can be supplied as family objects, zero-argument
#' constructors/factories, or names.
#' The name `"negative.binomial"` fits an estimated-size, log-link model with
#' `MASS::glm.nb()`. A `MASS::negative.binomial(theta, link)` object instead fits
#' a fixed-size model with `stats::glm()`. Common Gaussian, Gamma and inverse
#' Gaussian dispersions are not estimated for score preparation because their
#' common scale cancels from the empirically studentized statistic.
#'
#' @param formula A two-sided null-model formula.
#' @param data A data frame containing the model observations.
#' @param family A supported family object, zero-argument constructor/factory,
#'   or family name.
#' @param control A list of arguments accepted by `stats::glm.control()`.
#' @param tol Relative rank tolerance passed to `prepare_glm_score()`.
#' @return An internal `perglm_glm_score` object, including the fitted GLM in
#'   `$fit` and the prepared score residuals, weights and projection.
#' @noRd
fit_glm_null <- function(formula, data, family = stats::gaussian(),
                         control = stats::glm.control(), tol = 1e-7) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("formula must be a two-sided null-model formula.", call. = FALSE)
  }
  if (!is.data.frame(data) || nrow(data) < 2L) {
    stop("data must be a data frame containing at least two observations.",
         call. = FALSE)
  }
  if (!is.list(control)) {
    stop("control must be a list of arguments for stats::glm.control().",
         call. = FALSE)
  }
  control <- do.call(stats::glm.control, control)
  original_call <- match.call()
  estimated_nb <- is.character(family) && length(family) == 1L &&
    !is.na(family) && family == "negative.binomial"

  if (estimated_nb) {
    fit <- MASS::glm.nb(
      formula = formula, data = data, control = control,
      model = TRUE, x = TRUE, y = TRUE, na.action = stats::na.fail
    )
  } else {
    if (is.character(family)) {
      if (length(family) != 1L || is.na(family)) {
        stop("family must identify one supported GLM family.", call. = FALSE)
      }
      constructor <- switch(family,
        gaussian = stats::gaussian,
        binomial = stats::binomial,
        poisson = stats::poisson,
        Gamma = stats::Gamma,
        inverse.gaussian = stats::inverse.gaussian,
        NULL
      )
      if (is.null(constructor)) {
        stop("Unsupported GLM family name.", call. = FALSE)
      }
      family <- constructor()
    } else if (is.function(family)) {
      family <- family()
    }
    family <- standard_glm_family(family)$family
    fit <- stats::glm(
      formula = formula, data = data, family = family, control = control,
      model = TRUE, x = TRUE, y = TRUE, na.action = stats::na.fail,
      singular.ok = FALSE
    )
  }

  # Keep the caller's expressions rather than local wrapper variable names in
  # the ordinary GLM call. The retained model remains authoritative for scores.
  fit$call$formula <- original_call$formula
  fit$call$data <- original_call$data
  fit$call$control <- if (is.null(original_call$control)) {
    quote(stats::glm.control())
  } else {
    substitute(do.call(stats::glm.control, supplied_control),
               list(supplied_control = original_call$control))
  }
  if (!estimated_nb) {
    fit$call$family <- if (is.null(original_call$family)) {
      quote(stats::gaussian())
    } else original_call$family
  }

  prepared <- prepare_glm_score(fit, tol = tol)
  prepared$call <- original_call
  prepared
}
