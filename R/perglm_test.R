#' Permutation-calibrated score test for an exposure in a GLM
#'
#' Test for an exposure contribution to a conditional mean using an empirically
#' studentized GLM score and a fixed number of marginal exposure permutations.
#' The null model is fitted once; its working scores, weights and projection
#' factors are reused across permutations.
#'
#' @param formula A two-sided null-model formula containing the response and
#'   nuisance covariates only. Include an intercept and exclude the exposure.
#'   Fixed offsets can be included with `offset()`.
#' @param data A data frame containing the null-model observations.
#' @param exposure A finite real numeric vector with one value per row of
#'   `data`, in the same order. This argument is evaluated normally: use
#'   `exposure = dat$x`, rather than a bare column name or a character string.
#'   If named, names must exactly match the null-model row names.
#' @param family A supported GLM family object, zero-argument constructor or
#'   factory, or family name. The default is `stats::gaussian()`.
#' @param B A positive integer number of permutations; defaults to `999L`.
#' @param alternative One of `"two.sided"` (default), `"greater"`, or `"less"`.
#'   One-sided alternatives refer to the exposure coefficient on the link scale.
#' @param batch_size A positive integer maximum number of exposure permutations
#'   scored together. Smaller batches reduce working memory.
#' @param keep_statistics Whether to retain the permuted statistics in the
#'   result. Defaults to `FALSE`.
#' @param control A list of fitting controls accepted by [stats::glm.control()].
#' @param tol Relative tolerance for the nuisance-design rank check and
#'   residual exposure identifiability, strictly between zero and one.
#'
#' @details
#' The working alternative adds an exposure term to the linear predictor,
#' \eqn{g(\mu_i) = Z_i^T\beta + \gamma X_i + \mathrm{offset}_i}.
#' The null hypothesis is \eqn{\gamma = 0}. The alternatives `"greater"` and
#' `"less"` use the positive and negative signed-score tails, respectively.
#' With a decreasing link, a positive exposure coefficient corresponds to a
#' lower conditional mean.
#'
#' Supported working families are Gaussian, Bernoulli (`binomial`), Poisson,
#' Gamma, inverse Gaussian, and negative binomial, with their standard named
#' links. Binomial responses must be coded as zero or one. Grouped binomial
#' counts, prior weights, missing model observations, and rank-deficient or
#' intercept-free nuisance designs are not supported. The null fit must
#' converge to an interior solution and have positive residual degrees of
#' freedom. The exposure must have identifiable residual variation.
#' Bernoulli fits require both outcome values; count fits require at least one
#' positive outcome. A constant response with a constant offset gives zero
#' empirical score variance and is rejected.
#'
#' During null preparation, Bernoulli designs are checked once for complete or
#' quasi-complete separation using [detectseparation::detect_separation()].
#' The check uses the geometry of the free nuisance design, excluding fixed
#' offsets, and applies to every supported binomial link. Separated designs
#' cannot have an interior likelihood optimum and are rejected. The requested
#' working link and ordinary maximum-likelihood fitting procedure are preserved.
#'
#' Gaussian fits are checked once for residuals indistinguishable from fitting
#' roundoff. The numerical envelope follows the response, fitted means and
#' inverse-link derivative times absolute linear-predictor products. It scales
#' with measurement units and accounts for cancellation in the prediction.
#' This is a numerical-resolution diagnostic, not a statistical minimum-noise
#' threshold or a certificate for arbitrarily ill-conditioned fits.
#'
#' Supply `family = "negative.binomial"` to estimate the size parameter using
#' [MASS::glm.nb()] with a log link. Supply a
#' [MASS::negative.binomial()] family object to use a fixed size parameter.
#' Common dispersion factors cancel from empirical studentization.
#'
#' Each draw independently and uniformly permutes the raw exposure labels.
#' Repeated permutations and the identity permutation are allowed. Exposure
#' residuals and empirical score variance are recomputed for every draw while
#' the null-model quantities remain fixed. Inclusive tail comparisons give
#' \eqn{p = (1 + K)/(B + 1)}, where \eqn{K} is the exceedance count.
#' Two-sided comparisons use absolute scores, and the smallest attainable
#' p-value is \eqn{1/(B+1)}.
#'
#' Use `set.seed()` before calling the function for reproducibility. Successful
#' runs give identical results and final RNG states across batch sizes and
#' retention settings. If an observed or sampled score is undefined, including
#' zero empirical variance or numerical nuisance aliasing, the function stops
#' without a p-value. Draws are never discarded or redrawn.
#'
#' Conditional-mean inference follows the paper under a correct null mean model
#' and its regularity conditions; the working outcome variance may be
#' misspecified. The procedure does not test all forms of conditional dependence.
#'
#' @return An object of class `"htest"`, printed using the standard R method,
#'   with components:
#'   \describe{
#'     \item{statistic}{The signed empirically studentized score, named `T`.}
#'     \item{parameter}{The permutation budget, named `B`.}
#'     \item{p.value}{The permutation p-value.}
#'     \item{null.value}{The null exposure coefficient, zero.}
#'     \item{alternative}{The alternative hypothesis.}
#'     \item{method}{A description including the working family and link.}
#'     \item{data.name}{The exposure expression and null-model formula.}
#'     \item{B, exceedances}{The number of draws and their exceedance count.}
#'     \item{nobs}{The number of model observations.}
#'     \item{family, link, theta}{The working family name, link, and
#'       negative-binomial size (`NULL` for other families).}
#'     \item{null_diagnostics}{Null-fit convergence, iterations, and maximum
#'       normalized nuisance score, plus `separation_checked` and
#'       `gaussian_residual_check`. The latter contains `residual_norm`,
#'       `roundoff_bound`, and `relative_residual_norm` for Gaussian fits;
#'       otherwise it is `NULL`.}
#'     \item{permuted_statistics}{The sampled scores if `keep_statistics = TRUE`,
#'       otherwise `NULL`.}
#'   }
#'   No full model fit, coefficient estimate, or confidence interval is retained.
#'
#' @references
#' Ai, J. and Niu, Z. (2026). Conditional independence testing through marginal
#' permutation. \doi{10.48550/arXiv.2611.03255}.
#'
#' @examples
#' set.seed(42)
#' dat <- data.frame(z1 = rnorm(120), z2 = rnorm(120),
#'                   x = rbinom(120, 1, 0.3))
#' dat$y <- rpois(120, exp(0.3 + 0.2 * dat$z1 - 0.1 * dat$z2 + 0.2 * dat$x))
#'
#' set.seed(2026)
#' result <- perglm_test(y ~ z1 + z2, data = dat, exposure = dat$x,
#'                       family = stats::poisson(), B = 99L)
#' result
#' result$p.value
#'
#' @export
perglm_test <- function(formula, data, exposure, family = stats::gaussian(),
                         B = 999L,
                         alternative = c("two.sided", "greater", "less"),
                         batch_size = 64L, keep_statistics = FALSE,
                         control = stats::glm.control(), tol = 1e-7) {
  exposure_name <- deparse1(substitute(exposure))
  null <- fit_glm_null(formula = formula, data = data, family = family,
                       control = control, tol = tol)
  calibration <- calibrate_glm_permutation(
    null = null, exposure = exposure, B = B, alternative = alternative,
    batch_size = batch_size, keep_statistics = keep_statistics, tol = tol
  )
  structure(list(
    statistic = c(T = calibration$statistic),
    parameter = c(B = calibration$B),
    p.value = calibration$p.value,
    null.value = c("exposure coefficient" = 0),
    alternative = calibration$alternative,
    method = sprintf("Permutation-calibrated GLM score test (%s, %s link)",
                     null$family_name, null$family$link),
    data.name = paste(exposure_name, "with null model", deparse1(formula)),
    B = calibration$B,
    exceedances = calibration$exceedances,
    nobs = null$nobs,
    family = null$family_name,
    link = null$family$link,
    theta = null$theta,
    null_diagnostics = null$diagnostics,
    permuted_statistics = calibration$permuted_statistics
  ), class = "htest")
}
