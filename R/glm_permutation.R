#' Calibrate a GLM score with a fixed permutation budget
#'
#' Draw `B` independent, uniform permutations of the raw exposure labels.
#' Draws may repeat and may include the identity permutation. The null fit,
#' working score residuals, weights and projection factors stay fixed; each
#' exposure is residualized again and its score variance is recomputed.
#' The inclusive-tail Monte Carlo p-value is `(1 + exceedances) / (B + 1)`.
#'
#' Randomness comes only from `sample.int()`, so `set.seed()` controls the
#' result. On successful runs, changing `batch_size` or `keep_statistics` leaves
#' the draws and RNG state unchanged. The engine processes blocks of at most
#' `min(B, batch_size)` exposure columns; retaining scalar statistics is optional.
#'
#' Any nonidentifiable exposure, zero empirical variance or nonfinite statistic
#' stops the calculation without a p-value. Undefined draws are never discarded
#' or redrawn. This initial engine therefore requires well-defined scores on all
#' sampled relabelings; it does not define an extension for degenerate scores.
#'
#' @param null An object returned by `fit_glm_null()` or `prepare_glm_score()`.
#' @param exposure A finite real numeric vector in null-model row order.
#'   If named, names must exactly match the null-model rows.
#' @param B A positive integer permutation budget, at most `.Machine$integer.max`.
#' @param alternative One of `"two.sided"`, `"greater"`, or `"less"`.
#' @param batch_size A positive integer maximum number of exposure columns
#'   passed to the compiled batch kernel at once.
#' @param keep_statistics Whether to retain the `B` permuted statistics.
#' @param tol Relative exposure identifiability tolerance; see
#'   `glm_score_statistic()`.
#' @return An internal list with `statistic`, `p.value`, `alternative`, `B`,
#'   `exceedances`, and `permuted_statistics` (`NULL` unless retained).
#' @noRd
calibrate_glm_permutation <- function(null, exposure, B = 999L,
                                      alternative = c("two.sided", "greater", "less"),
                                      batch_size = 64L, keep_statistics = FALSE,
                                      tol = 1e-7) {
  validate_glm_score_inputs(null, exposure, tol)
  B <- permutation_positive_integer(B, "B")
  batch_size <- permutation_positive_integer(batch_size, "batch_size")
  alternative <- match.arg(alternative)
  if (!is.logical(keep_statistics) || length(keep_statistics) != 1L ||
      !is.null(dim(keep_statistics)) || is.na(keep_statistics)) {
    stop("keep_statistics must be TRUE or FALSE.", call. = FALSE)
  }
  observed <- glm_score_statistic(null, exposure, tol = tol)
  n <- null$nobs
  # Labels carry their values into fixed design rows. Input names have already
  # been checked and must not be interpreted as new row identities after a draw.
  exposure <- unname(exposure)
  permuted_statistics <- if (keep_statistics) numeric(B) else NULL
  exposures <- matrix(NA_real_, nrow = n, ncol = min(batch_size, B))
  exceedances <- 0
  completed <- 0L
  while (completed < B) {
    columns <- min(batch_size, B - completed)
    if (ncol(exposures) != columns) {
      exposures <- matrix(NA_real_, nrow = n, ncol = columns)
    }
    for (column in seq_len(columns)) {
      exposures[, column] <- exposure[sample.int(n, n, replace = FALSE)]
    }
    statistics <- tryCatch(
      residualized_score_batch_cpp(null$projection, exposures,
        null$score_residuals, exposure_tol = tol),
      error = function(error) {
        stop(sprintf("Permutations %d-%d: %s", completed + 1L,
                     completed + columns, conditionMessage(error)), call. = FALSE)
      }
    )
    if (any(!is.finite(statistics))) {
      bad_draw <- completed + which(!is.finite(statistics))[1L]
      stop(sprintf("Permutation %d: the GLM score statistic must be finite.", bad_draw),
           call. = FALSE)
    }
    exceedances <- exceedances + sum(permutation_extremes(observed, statistics,
                                                          alternative))
    if (keep_statistics) {
      permuted_statistics[completed + seq_len(columns)] <- statistics
    }
    completed <- completed + columns
  }
  list(statistic = observed, p.value = (1 + exceedances) / (B + 1),
       alternative = alternative, B = B, exceedances = exceedances,
       permuted_statistics = permuted_statistics)
}

permutation_positive_integer <- function(value, name) {
  if (!glm_numeric_vector(value, 1L) || value < 1 ||
      value > .Machine$integer.max || value != floor(value)) {
    stop(sprintf("%s must be a positive integer no larger than .Machine$integer.max.",
                 name), call. = FALSE)
  }
  as.integer(value)
}
