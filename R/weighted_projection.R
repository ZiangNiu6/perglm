#' Prepare a weighted projection
#'
#' Cache a QR decomposition and its thin factors for repeated weighted
#' regressions on the same design in R or compiled kernels.
#' This numerical primitive does not require an intercept.
#'
#' @param z A finite real numeric design matrix with full column rank.
#' @param weights Finite, strictly positive weights, one per row of `z`.
#' @param tol Relative tolerance used by the QR rank check.
#' @return An internal `perglm_projection` object.
#' @noRd
prepare_projection <- function(z, weights, tol = 1e-7) {
  if (!is.matrix(z) || !is.numeric(z) || is.complex(z) ||
      nrow(z) < 1L || ncol(z) < 1L || any(!is.finite(z))) {
    stop("z must be a nonempty, finite real numeric matrix.", call. = FALSE)
  }
  if (nrow(z) < ncol(z)) {
    stop("z must have at least as many rows as columns.", call. = FALSE)
  }
  if (!is.numeric(weights) || is.complex(weights) ||
      !is.null(dim(weights)) || length(weights) != nrow(z) ||
      any(!is.finite(weights)) || any(weights <= 0)) {
    stop("weights must be a finite, strictly positive vector with one value per row of z.",
         call. = FALSE)
  }
  if (!is.numeric(tol) || is.complex(tol) || length(tol) != 1L ||
      !is.finite(tol) || tol <= 0 || tol >= 1) {
    stop("tol must be a finite number strictly between zero and one.",
         call. = FALSE)
  }

  # A common weight scale leaves the projection unchanged. Taking square roots
  # before normalizing avoids prematurely underflowing very small weights.
  sqrt_weights <- sqrt(weights)
  sqrt_weights <- sqrt_weights / max(sqrt_weights)
  weighted_qr <- qr(sqrt_weights * z, tol = tol, LAPACK = FALSE)
  if (weighted_qr$rank != ncol(z)) {
    stop("z must have full column rank under the supplied weights and tolerance.",
         call. = FALSE)
  }

  structure(
    list(z = z, sqrt_weights = sqrt_weights, qr = weighted_qr,
         q = qr.Q(weighted_qr, complete = FALSE), r = qr.R(weighted_qr)),
    class = "perglm_projection"
  )
}

#' Residualize an exposure using a prepared projection
#'
#' @param projection An object returned by `prepare_projection()`.
#' @param exposure A finite real numeric vector, aligned with the design rows.
#' @return The unweighted residual vector from weighted least squares.
#' @noRd
residualize_exposure <- function(projection, exposure) {
  if (!inherits(projection, "perglm_projection")) {
    stop("projection must be an object returned by prepare_projection().",
         call. = FALSE)
  }
  if (!is.numeric(exposure) || is.complex(exposure) ||
      !is.null(dim(exposure)) || length(exposure) != nrow(projection$z) ||
      any(!is.finite(exposure))) {
    stop("exposure must be a finite real numeric vector with one value per design row.",
         call. = FALSE)
  }

  coefficients <- qr.coef(
    projection$qr,
    projection$sqrt_weights * exposure
  )
  # Computing residuals on the original scale avoids division by tiny weights.
  residuals <- exposure - drop(projection$z %*% coefficients)
  if (any(!is.finite(residuals))) {
    stop("Exposure residualization produced nonfinite values.", call. = FALSE)
  }
  residuals
}
