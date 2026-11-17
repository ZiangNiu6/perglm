#' Empirically studentize score contributions
#'
#' Compute `sum(r) / sqrt(sum((r - mean(r))^2))`, corresponding to the
#' empirical variance with divisor `length(r)`. Constant contributions have
#' zero variance and cause an error. No near-zero variance threshold is imposed.
#'
#' @param contributions A finite real numeric vector of at least two scores.
#' @return A signed scalar score statistic.
#' @noRd
studentize_score <- function(contributions) {
  if (!is.numeric(contributions) || is.complex(contributions) ||
      !is.null(dim(contributions)) || length(contributions) < 2L ||
      any(!is.finite(contributions))) {
    stop("contributions must be a finite real numeric vector of length at least two.",
         call. = FALSE)
  }

  # Scale the numerator to avoid overflow. For the denominator, subtract an
  # original-scale anchor first so tiny differences between near-constant
  # contributions survive normalization.
  scale <- max(abs(contributions))
  if (scale == 0) {
    stop("Score contributions must have positive empirical variance.",
         call. = FALSE)
  }
  scaled <- contributions / scale
  differences <- contributions - contributions[1L]
  if (all(is.finite(differences))) {
    difference_scale <- max(abs(differences))
    if (difference_scale == 0) {
      stop("Score contributions must have positive empirical variance.",
           call. = FALSE)
    }
    scaled_differences <- differences / difference_scale
    centered <- scaled_differences - mean(scaled_differences)
    relative_scale <- difference_scale / scale
  } else {
    # Opposite-sign extremes can overflow anchor subtraction; their scaled
    # contributions have a substantial range and can be centered directly.
    centered <- scaled - mean(scaled)
    relative_scale <- 1
  }
  sum_squares <- sum(centered * centered)
  if (sum_squares == 0) {
    stop("Score contributions must have positive empirical variance.",
         call. = FALSE)
  }

  sum(scaled) / (relative_scale * sqrt(sum_squares))
}
