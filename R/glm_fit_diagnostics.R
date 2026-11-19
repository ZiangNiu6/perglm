#' Diagnose numerically unresolved Gaussian residual variation
#'
#' Called after the retained Gaussian fit and working link derivative have been
#' validated. The comparison is relative to response and prediction arithmetic,
#' rather than to an absolute minimum residual or a statistical noise threshold.
#'
#' For unit roundoff `u`, use `gamma = k*u/(1-k*u)`, with `k = 4*(p+2)`.
#' The normwise envelope is `gamma` times the sum of the norms of `y`, `mu`,
#' and `abs(derivative) * (abs(z) %*% abs(coefficients) + abs(offset))`.
#' The last term transfers linear-predictor evaluation uncertainty through the
#' inverse link. Its absolute products account for cancellation and are
#' unchanged by rescaling individual design columns and their coefficients.
#'
#' This conservative envelope represents a small number of floating-point
#' operations per predictor term. It is a numerical-resolution diagnostic, not
#' a certified global error bound for QR or nonlinear GLM iteration. In
#' particular, arbitrarily ill-conditioned fits are not certified by passing
#' this check. Fits whose residuals are indistinguishable from roundoff are
#' unsupported; very small but numerically resolved residual noise is allowed.
#'
#' @param y,mu Validated finite response and final fitted-mean vectors.
#' @param z Validated finite nuisance-design matrix.
#' @param coefficients Validated finite coefficients in design-column order.
#' @param offset Validated finite offset vector, including zeros when absent.
#' @param derivative The validated finite nonzero inverse-link derivative.
#' @return A list with residual norm, roundoff bound, and residual norm relative
#'   to the sum of response and fitted-mean norms. Norms use the response units.
#' @noRd
gaussian_fit_diagnostics <- function(y, mu, z, coefficients, offset, derivative) {
  response_scale <- max(abs(y), abs(mu))
  if (response_scale == 0) {
    stop("The Gaussian null fit has no numerically resolved residual variation; residuals are indistinguishable from fitting roundoff.",
         call. = FALSE)
  }

  scaled_y <- y / response_scale
  scaled_mu <- mu / response_scale
  residuals <- y - mu
  scaled_residuals <- if (all(is.finite(residuals))) {
    residuals / response_scale
  } else {
    # Subtracting opposite-sign extremes may overflow on the original scale.
    scaled_y - scaled_mu
  }
  residual_norm <- gaussian_diagnostic_norm(scaled_residuals)
  response_norm <- gaussian_diagnostic_norm(scaled_y) +
    gaussian_diagnostic_norm(scaled_mu)

  predictor_magnitude <- drop(abs(z) %*% abs(coefficients)) + abs(offset)
  mapped_magnitude <- numeric(length(y))
  positive <- predictor_magnitude > 0
  # Logarithms are used only to scale the error envelope, never the residuals
  # or score. They avoid intermediate products overflowing or underflowing.
  mapped_magnitude[positive] <- exp(log(abs(derivative[positive])) +
    log(predictor_magnitude[positive]) - log(response_scale))

  unit_roundoff <- .Machine$double.eps / 2
  operations <- 4 * (ncol(z) + 2)
  gamma <- operations * unit_roundoff / (1 - operations * unit_roundoff)
  roundoff_bound <- gamma *
    (response_norm + gaussian_diagnostic_norm(mapped_magnitude))
  if (residual_norm <= roundoff_bound) {
    stop("The Gaussian null fit has no numerically resolved residual variation; residuals are indistinguishable from fitting roundoff.",
         call. = FALSE)
  }

  list(residual_norm = response_scale * residual_norm,
       roundoff_bound = response_scale * roundoff_bound,
       relative_residual_norm = residual_norm / response_norm)
}

#' @noRd
gaussian_diagnostic_norm <- function(value) {
  scale <- max(abs(value))
  if (scale == 0 || is.infinite(scale)) {
    return(scale)
  }
  scale * sqrt(sum((value / scale)^2))
}
