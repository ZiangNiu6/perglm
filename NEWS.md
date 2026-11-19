# perglm 0.0.0.9000

* Introduce `perglm_test()` for permutation-calibrated, empirically studentized
  GLM score tests, returning a standard R `htest` object.
* Support Gaussian, Bernoulli, Poisson, Gamma, inverse Gaussian, and negative
  binomial working families, with standard named links and fixed offsets.
  Negative binomial size can be fixed or estimated under a log-link null model.
* Provide fixed-budget marginal exposure permutations, two-sided and one-sided
  alternatives, inclusive ties, and the plus-one p-value correction. Support
  reproducibility with `set.seed()` and optional retention of permuted scores.
* Fit the nuisance-only null model once and reuse score quantities and QR
  factors. Use compiled C++ kernels and bounded batches for permutation scoring.
* Require complete observations, unit prior weights, a full-rank nuisance design
  with an intercept, and identifiable exposure variation. Check null-model
  convergence and empirical score variance. Reject separated Bernoulli fits and
  Gaussian fits with residuals indistinguishable from numerical roundoff.
  Undefined sampled scores stop the test without discarding or redrawing
  permutations.
* Add an executable getting-started vignette, a package website, installation
  instructions, and citations for the method and software.
* Add GitHub Actions checks on Linux, macOS, and Windows, including vignettes
  and the minimum supported R version, 4.1.0.
* Require testthat 3.1.7 or later for the namespace mocking used by the tests.
