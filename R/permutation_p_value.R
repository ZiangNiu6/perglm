#' Calculate a fixed-budget Monte Carlo permutation p-value
#'
#' Count inclusive ties in the requested tail and include the observed statistic
#' through the plus-one correction. Two-sided comparisons use absolute values.
#' This helper performs no sampling and does not consume random numbers.
#'
#' @param observed A finite real numeric scalar statistic.
#' @param permuted A nonempty finite real numeric vector of permuted statistics.
#' @param alternative The alternative hypothesis.
#' @return An unnamed scalar p-value in `(0, 1]`.
#' @noRd
permutation_p_value <- function(observed, permuted,
                                alternative = c("two.sided", "greater", "less")) {
  if (!is.numeric(observed) || is.complex(observed) ||
      !is.null(dim(observed)) || length(observed) != 1L ||
      !is.finite(observed)) {
    stop("observed must be a finite real numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(permuted) || is.complex(permuted) ||
      !is.null(dim(permuted)) || length(permuted) < 1L ||
      any(!is.finite(permuted))) {
    stop("permuted must be a nonempty finite real numeric vector.",
         call. = FALSE)
  }
  alternative <- match.arg(alternative)
  extremes <- permutation_extremes(observed, permuted, alternative)
  (1 + sum(extremes)) / (length(permuted) + 1)
}

#' Compare already validated observed and permutation statistics
#'
#' @param observed A finite real numeric scalar.
#' @param permuted A finite real numeric vector.
#' @param alternative One of `"two.sided"`, `"greater"`, or `"less"`.
#' @return A logical vector indicating inclusive tail comparisons.
#' @noRd
permutation_extremes <- function(observed, permuted, alternative) {
  switch(alternative,
    two.sided = abs(permuted) >= abs(observed),
    greater = permuted >= observed,
    less = permuted <= observed
  )
}
