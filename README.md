<!-- README.md is generated from README.Rmd. Please edit that file -->

# perglm <a href="https://github.com/ZiangNiu6/perglm"><img src="man/figures/logo.png" align="right" height="139" alt="perglm hex logo" /></a>

<!-- badges: start -->
<!-- badges: end -->

**perglm** is being developed to provide permutation-calibrated, empirically
studentized score tests for generalized linear models.

## Development status

Internal helpers fit null GLMs and prepare their working score residuals,
weights, and cached projections. Supported families are Gaussian, Bernoulli,
Poisson, Gamma, inverse Gaussian, and negative binomial, including fixed offsets.
Model observations must be complete, with unit prior weights and a full-rank
nuisance design containing an intercept.

Rcpp kernels support scalar and batch studentization, as well as combined
exposure residualization and score calculation, alongside an R reference
implementation. The public GLM score-test interface and permutation engines are
under development.

The package layout includes R code in `R/`, native code reserved for `src/`,
documentation in `man/`, tests in `tests/`, and pkgdown configuration for a future
website in `docs/`.

## Reference

Ai, J. and Niu, Z. (2026). *Conditional independence testing through marginal
permutation*. [arXiv:2611.03255](https://arxiv.org/abs/2611.03255).

## License

GPL version 3 or later. See [LICENSE.md](LICENSE.md).
