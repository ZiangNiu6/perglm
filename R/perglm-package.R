#' Permutation Score Tests for Generalized Linear Models
#'
#' Package infrastructure for permutation-calibrated, empirically studentized
#' score tests in generalized linear models. Internal helpers fit null models,
#' prepare working scores and weights, and provide weighted exposure
#' residualization and empirical studentization. Numerical calculations have
#' R reference implementations and compiled scalar and batch kernels. Internal
#' calibration helpers provide fixed-budget permutation tests with inclusive
#' tails and a plus-one correction. The public testing interface is under
#' development.
#'
#' @references
#' Ai, J. and Niu, Z. (2026). Conditional independence testing through marginal
#' permutation. \doi{10.48550/arXiv.2611.03255}.
#'
#' @useDynLib perglm, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @keywords internal
"_PACKAGE"
