.separation_fixture <- function(kind = c("overlap", "complete", "quasi")) {
  kind <- match.arg(kind)
  i <- seq_len(40)
  z <- rep(c(-1, 1), each = 20)
  y <- switch(kind,
    overlap = rep(0:1, 20),
    complete = as.integer(z > 0),
    quasi = c(rep(0:1, 10), rep(1, 20))
  )
  data.frame(y = y, z = z,
    x = sin(i * 0.71) + cos(i * 0.17),
    off = 0.1 * sin(i * 0.37),
    row.names = paste0("case", i))
}

test_that("complete and quasi separation fail before permutation sampling", {
  for (link in c("logit", "probit", "cloglog", "cauchit", "log", "identity")) {
    for (kind in c("complete", "quasi")) {
      dat <- .separation_fixture(kind)
      set.seed(17)
      before <- .Random.seed
      # Bounded-domain links can fail during fitting before the separation
      # guard. Either route must prevent inference from a noninterior fit.
      expect_error(suppressWarnings(perglm_test(y ~ z, dat, dat$x,
        family = stats::binomial(link), B = 9L,
        control = list(maxit = 200L))),
        "separation|finite|interior|converg|valid")
      expect_identical(.Random.seed, before, info = paste(link, kind))
    }
  }
})

test_that("separation detection respects nuisance bases and factor reference levels", {
  for (kind in c("complete", "quasi")) {
    dat <- .separation_fixture(kind)
    dat$group <- factor(ifelse(dat$z < 0, "left", "right"),
                         levels = c("left", "right"))
    for (formula in list(y ~ z, y ~ I(17 + 0.125 * z), y ~ group)) {
      expect_error(perglm_test(formula, dat, dat$x, family = stats::binomial(),
        B = 1L, control = list(maxit = 200L)), "separation|finite")
    }
    dat$group <- stats::relevel(dat$group, ref = "right")
    expect_error(perglm_test(y ~ group, dat, dat$x, family = stats::binomial(),
      B = 1L, control = list(maxit = 200L)), "separation|finite")
  }
})

test_that("full-rank overlapping responses remain valid across working links", {
  dat <- .separation_fixture("overlap")
  set.seed(819)
  for (draw in seq_len(5)) sample.int(nrow(dat), nrow(dat), replace = FALSE)
  expected_seed <- .Random.seed
  for (link in c("logit", "probit", "cloglog", "cauchit", "log", "identity")) {
    set.seed(819)
    result <- perglm_test(y ~ z, dat, dat$x, family = stats::binomial(link),
      B = 5L, control = list(maxit = 200L))
    expect_s3_class(result, "htest")
    expect_identical(result$link, link)
    expect_true(is.finite(result$statistic))
    expect_true(result$p.value >= 1 / 6 && result$p.value <= 1)
    expect_identical(.Random.seed, expected_seed, info = link)
  }
})

test_that("known offsets can classify outcomes without free-coefficient separation", {
  dat <- .separation_fixture("overlap")
  dat$off <- ifelse(dat$y == 1, 2, -2)
  for (link in c("logit", "probit", "cloglog", "cauchit")) {
    fit <- stats::glm(y ~ 1 + offset(off), data = dat,
      family = stats::binomial(link), control = stats::glm.control(maxit = 200L))
    expect_true(all((2 * dat$y - 1) * fit$linear.predictors > 0))
    set.seed(17)
    result <- perglm_test(y ~ 1 + offset(off), dat, dat$x,
      family = stats::binomial(link), B = 9L,
      control = list(maxit = 200L))
    expect_s3_class(result, "htest")
    expect_true(is.finite(result$p.value))
  }
})

test_that("finite varying offsets do not remove an actual separating direction", {
  for (link in c("logit", "probit", "cloglog", "cauchit")) {
    dat <- .separation_fixture("complete")
    set.seed(17)
    before <- .Random.seed
    expect_error(suppressWarnings(perglm_test(y ~ z + offset(off), dat, dat$x,
      family = stats::binomial(link), B = 9L,
      control = list(maxit = 200L))), "separation|finite|interior|converg|valid")
    expect_identical(.Random.seed, before, info = link)
  }
})

test_that("public outcomes agree with an external geometric separation detector", {
  skip_if_not_installed("detectseparation")
  for (kind in c("overlap", "complete", "quasi")) {
    dat <- .separation_fixture(kind)
    design <- stats::model.matrix(y ~ z, data = dat)
    detector <- detectseparation::detect_separation(
      x = design, y = dat$y, offset = dat$off,
      family = stats::binomial("logit"), intercept = TRUE
    )
    separated <- any(is.infinite(detector$coefficients))
    expect_identical(separated, kind != "overlap")
    set.seed(17)
    result <- tryCatch(perglm_test(y ~ z + offset(off), dat, dat$x,
      family = stats::binomial(), B = 3L,
      control = list(maxit = 200L)), error = identity)
    expect_identical(inherits(result, "error"), separated)
  }
})

test_that("the expensive detector is called once per null fit rather than per draw", {
  skip_if_not_installed("detectseparation")
  detector <- getExportedValue("detectseparation", "detect_separation")
  calls <- 0L
  arguments <- NULL
  counted_detector <- function(...) {
    calls <<- calls + 1L
    arguments <<- list(...)
    detector(...)
  }
  testthat::local_mocked_bindings(detect_separation = counted_detector,
                                  .package = "detectseparation")
  dat <- .separation_fixture("overlap")
  dat$y <- factor(ifelse(dat$y == 1, "yes", "no"), levels = c("no", "yes"))
  for (link in c("logit", "log", "identity")) {
    for (B in c(1L, 9L)) {
      calls <- 0L
      arguments <- NULL
      set.seed(17)
      result <- perglm_test(y ~ z + offset(off), dat, dat$x,
        family = stats::binomial(link), B = B, batch_size = 2L,
        control = list(maxit = 200L))
      expect_identical(calls, 1L, info = paste(link, B))
      expect_identical(arguments$family$link, "logit")
      expect_equal(unname(arguments$y), as.integer(dat$y) - 1L)
      expect_identical(nrow(arguments$x), nrow(dat))
      expect_identical(result$link, link)
    }
  }
})
