#' Permutation Score Tests for Generalized Linear Models
#'
#' Use [perglm_test()] for permutation-calibrated, empirically studentized
#' score tests in generalized linear models. A nuisance-only null model is
#' fitted once and its score quantities are reused across fixed-budget exposure
#' permutations. Compiled scalar and batch kernels provide weighted exposure
#' residualization and empirical studentization, with inclusive tails and a
#' plus-one correction.
#'
#' @references
#' Ai, J. and Niu, Z. (2026). Conditional independence testing through marginal
#' permutation. \doi{10.48550/arXiv.2611.03255}.
#'
#' @useDynLib perglm, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @keywords internal
"_PACKAGE"
