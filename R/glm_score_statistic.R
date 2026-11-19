#' Calculate an empirically studentized GLM score
#'
#' Residualize the exposure using the cached null-model projection, multiply
#' by the fixed working score residuals, and calculate
#' `sum(r) / sqrt(sum((r - mean(r))^2))` in the compiled kernel.
#'
#' The exposure is centered and scaled internally by a positive affine
#' transformation, which leaves this statistic unchanged with an intercept.
#' This avoids sensitivity to exposure units and large common offsets.
#' An exposure is rejected when its weighted residual norm is at most `tol`
#' times its weighted norm after normalization. Zero empirical score variance
#' and nonfinite statistics also cause errors.
#'
#' @param null An object returned by `fit_glm_null()` or `prepare_glm_score()`.
#' @param exposure A finite real numeric vector in null-model row order.
#'   If named, its names must exactly match the retained row names.
#' @param tol Relative tolerance for residual exposure identifiability,
#'   strictly between zero and one.
#' @return A signed scalar score statistic. No null model is refitted.
#' @noRd
glm_score_statistic <- function(null, exposure, tol = 1e-7) {
  validate_glm_score_inputs(null, exposure, tol)
  statistic <- residualized_score_cpp(
    null$projection, exposure, null$score_residuals, exposure_tol = tol
  )
  if (!is.finite(statistic)) {
    stop("The GLM score statistic must be finite.", call. = FALSE)
  }
  statistic
}

# Shared by scalar scoring and permutation calibration; the native kernel
# separately validates the cached factors before using them.
validate_glm_score_inputs <- function(null, exposure, tol) {
  invalid_null <- function() {
    stop("null must be a valid object returned by fit_glm_null() or prepare_glm_score().",
         call. = FALSE)
  }
  if (!inherits(null, "perglm_glm_score") || !is.list(null) ||
      !is.list(null$projection) ||
      !inherits(null$projection, "perglm_projection") ||
      !glm_numeric_vector(null$nobs, 1L) || null$nobs < 2 ||
      null$nobs != floor(null$nobs) || !is.matrix(null$projection$z)) {
    invalid_null()
  }
  n <- null$nobs
  if (nrow(null$projection$z) != n ||
      !glm_numeric_vector(null$score_residuals, n) ||
      !is.character(null$row_names) || length(null$row_names) != n ||
      anyNA(null$row_names) || anyDuplicated(null$row_names) ||
      !identical(rownames(null$projection$z), null$row_names) ||
      (!is.null(names(null$score_residuals)) &&
       !identical(names(null$score_residuals), null$row_names))) {
    invalid_null()
  }
  if (!glm_numeric_vector(exposure, n)) {
    stop("exposure must be a finite real numeric vector with one value per null-model row.",
         call. = FALSE)
  }
  if (!is.null(names(exposure)) &&
      !identical(names(exposure), null$row_names)) {
    stop("Exposure names must exactly match the null-model row order.", call. = FALSE)
  }
  if (!glm_numeric_vector(tol, 1L) || tol <= 0 || tol >= 1) {
    stop("tol must be a finite number strictly between zero and one.", call. = FALSE)
  }
  invisible(NULL)
}
