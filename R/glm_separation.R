#' Check a Bernoulli nuisance design for separation
#'
#' Called once during null preparation, before any exposure permutations.
#' A separating direction in the free design strictly increases the Bernoulli
#' likelihood locally at every interior fit for any supported (increasing)
#' binomial link. Thus geometric separation precludes an interior solution,
#' including for log/identity links whose likelihood maximum can be finite on
#' the boundary. Fixed offsets are not free directions.
#'
#' The detector's logit family requests this geometry check only; the actual
#' working family and its maximum-likelihood fit are not changed. Positive
#' column normalization preserves separation and avoids dependence on units.
#'
#' @param z Validated, finite, full-rank nuisance design including an intercept.
#' @param y Validated Bernoulli response with both outcome classes present.
#' @param offset Validated fixed offsets, excluded from separating directions.
#' @param family The validated working binomial family.
#' @return Invisibly `TRUE` if the detector establishes no separation.
#' @noRd
check_glm_separation <- function(z, y, offset, family) {
  scales <- apply(abs(z), 2L, max)
  if (any(!is.finite(scales)) || any(scales <= 0)) {
    stop("The Bernoulli nuisance design must have full column rank.", call. = FALSE)
  }
  design <- sweep(z, 2L, scales, "/")
  detected <- tryCatch(
    withCallingHandlers(
      detectseparation::detect_separation(
        x = design, y = y, family = stats::binomial("logit"),
        intercept = TRUE, singular.ok = FALSE
      ),
      warning = function(warning) stop(conditionMessage(warning), call. = FALSE)
    ),
    error = function(error) {
      stop("Could not validate Bernoulli separation: ", conditionMessage(error),
           call. = FALSE)
    }
  )
  if (!is.list(detected) || !is.logical(detected$outcome) ||
      length(detected$outcome) != 1L || is.na(detected$outcome)) {
    stop("Could not establish whether the Bernoulli null model has separation.",
         call. = FALSE)
  }
  if (detected$outcome) {
    stop("The Bernoulli null model has complete or quasi-complete separation; a finite interior fit is required.",
         call. = FALSE)
  }
  invisible(TRUE)
}
