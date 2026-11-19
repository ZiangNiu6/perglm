<!-- README.md is generated from README.Rmd. Please edit that file -->

# perglm <a href="https://github.com/ZiangNiu6/perglm"><img src="man/figures/logo.png" align="right" height="139" alt="perglm hex logo" /></a>

<!-- badges: start -->
<!-- badges: end -->

**perglm** provides permutation-calibrated, empirically studentized score tests
for an exposure in generalized linear models.

## Quick start

Supply the nuisance-only null model as the formula, and pass the exposure
separately as a numeric vector. The result is a standard R `htest` object.

```r
library(perglm)

set.seed(42)
dat <- data.frame(z1 = rnorm(120), z2 = rnorm(120),
                  x = rbinom(120, 1, 0.3))
dat$y <- rpois(120, exp(0.3 + 0.2 * dat$z1 - 0.1 * dat$z2 + 0.2 * dat$x))

set.seed(2026)
result <- perglm_test(
  y ~ z1 + z2,
  data = dat,
  exposure = dat$x,
  family = poisson(),
  B = 999L
)
result
result$p.value
```

The default alternative is `"two.sided"`; `"greater"` and `"less"` refer to
the exposure coefficient on the link scale. Set `keep_statistics = TRUE` to
retain the sampled scores in `result$permuted_statistics`. See `?perglm_test`
for supported families, assumptions, and fitting controls.

The getting-started vignette walks through this example and its output:

```r
vignette("getting-started", package = "perglm")
```

Read the [rendered getting-started tutorial](https://ziangniu6.github.io/perglm/articles/getting-started.html)
on the package website. Its [R Markdown source](https://github.com/ZiangNiu6/perglm/blob/main/vignettes/getting-started.Rmd)
is in the source repository.

## Reference

Ai, J. and Niu, Z. (2026). *Conditional independence testing through marginal
permutation*. [arXiv:2611.03255](https://arxiv.org/abs/2611.03255).

## License

GPL version 3 or later. See [LICENSE.md](LICENSE.md).
