test_that("the fused score matches the composed R reference", {
  z <- cbind(1, c(-2, -1, 0, 1, 2, 3), c(0, 1, 0, -1, 1, 2))
  weights <- c(0.2, 1, 4, 0.5, 2, 3)
  projection <- prepare_projection(z, weights)
  score_residuals <- c(-2, 0.5, 1, -1.5, 2, 0.75)
  for (exposure in list(c(0L, 1L, 0L, 0L, 1L, 0L),
                        c(2, -3, 1, 5, 0, 4))) {
    reference <- studentize_score(
      residualize_exposure(projection, exposure) * score_residuals
    )
    fused <- residualized_score_cpp(projection, exposure, score_residuals)
    expect_lte(abs(fused - reference), 1e-11 * max(1, abs(reference)))
  }
})

test_that("fused scores correctly unmap pivoted coefficients", {
  z <- cbind(intercept = 1, large = 10 * c(-2, -1, 0, 1, 2, 3),
             small = c(0, 0.5, 0, -0.5, 0.5, 1))
  projection <- prepare_projection(z, c(0.2, 1, 4, 0.5, 2, 3))
  # The regular full-rank LINPACK cache need not pivot. Use equivalent LAPACK
  # factors to exercise a genuine nonidentity pivot without changing the model.
  projection$qr <- qr(projection$sqrt_weights * z, LAPACK = TRUE)
  projection$q <- qr.Q(projection$qr, complete = FALSE)
  projection$r <- qr.R(projection$qr)
  expect_false(identical(projection$qr$pivot, seq_len(ncol(z))))
  exposure <- c(2, -3, 1, 5, 0, 4)
  score_residuals <- c(-2, 0.5, 1, -1.5, 2, 0.75)
  reference <- studentize_score(
    residualize_exposure(projection, exposure) * score_residuals
  )
  expect_lte(abs(residualized_score_cpp(projection, exposure, score_residuals) -
                   reference), 1e-11 * max(1, abs(reference)))
})

test_that("batch fused scores match individual calls and retain names", {
  set.seed(627)
  z <- cbind(1, rnorm(40), rnorm(40))
  projection <- prepare_projection(z, runif(40, 0.1, 3))
  exposures <- cbind(binary = rbinom(40, 1, 0.3),
                     normal = rnorm(40), scaled = 100 * rnorm(40))
  score_residuals <- rnorm(40)
  original_exposures <- exposures
  original_residuals <- score_residuals
  original_projection <- projection
  reference <- apply(exposures, 2L, function(exposure) {
    studentize_score(residualize_exposure(projection, exposure) * score_residuals)
  })
  individual <- vapply(seq_len(ncol(exposures)), function(column) {
    residualized_score_cpp(projection, exposures[, column], score_residuals)
  }, numeric(1))
  batch <- residualized_score_batch_cpp(projection, exposures, score_residuals)
  expect_equal(batch, reference, tolerance = 1e-11)
  expect_identical(unname(batch), individual)
  expect_identical(exposures, original_exposures)
  expect_identical(score_residuals, original_residuals)
  expect_identical(projection, original_projection)
  expect_identical(
    residualized_score_batch_cpp(projection, matrix(numeric(), 40, 0),
                                  score_residuals), numeric()
  )
  expect_equal(
    residualized_score_batch_cpp(projection, exposures[, 1L, drop = FALSE],
                                  score_residuals),
    c(binary = individual[1L]), tolerance = 1e-14
  )
})

test_that("integer designs and exposure blocks are supported", {
  z <- cbind(1L, -2L:3L)
  projection <- prepare_projection(z, rep(1L, 6))
  exposures <- cbind(c(0L, 1L, 0L, 0L, 1L, 0L), c(2L, -3L, 1L, 5L, 0L, 4L))
  score_residuals <- c(-2L, 1L, 1L, -1L, 2L, 1L)
  reference <- apply(exposures, 2L, function(exposure) {
    studentize_score(residualize_exposure(projection, exposure) * score_residuals)
  })
  expect_equal(residualized_score_batch_cpp(projection, exposures, score_residuals),
               reference, tolerance = 1e-11)
})

test_that("invalid fused score inputs are rejected", {
  projection <- prepare_projection(cbind(1, c(-2, -1, 0, 1, 2, 3)), rep(1, 6))
  exposure <- c(0, 1, 0, 0, 1, 0)
  score_residuals <- c(-2, 1, 1, -1, 2, 1)
  expect_error(residualized_score_cpp(NULL, exposure, score_residuals),
               "projection must")
  for (invalid in list(1:5, c(NA_integer_, 1:5), c(Inf, 1:5),
                       rep(TRUE, 6), rep(1i, 6), matrix(1:6, 6),
                       factor(letters[1:6]))) {
    expect_error(residualized_score_cpp(projection, invalid, score_residuals),
                 "exposure must")
    expect_error(residualized_score_cpp(projection, exposure, invalid),
                 "score_residuals must")
  }
  for (invalid in list(exposure, matrix(1:5, 5), matrix(TRUE, 6, 2),
                       matrix(1i, 6, 2))) {
    expect_error(residualized_score_batch_cpp(projection, invalid, score_residuals),
                 "exposures must")
  }
  expect_error(residualized_score_batch_cpp(projection,
                 cbind(exposure, c(NA_integer_, 1:5)), score_residuals),
               "Column 2:.*exposures must")
  expect_error(residualized_score_cpp(projection, exposure, rep(0, 6)),
               "positive empirical variance")
  expect_error(residualized_score_batch_cpp(projection, cbind(exposure, exposure),
                                            rep(0, 6)),
               "Column 1:.*positive empirical variance")
})

test_that("incompatible or corrupted projection caches are rejected", {
  projection <- prepare_projection(cbind(1, c(-2, -1, 0, 1, 2, 3)), rep(1, 6))
  exposure <- c(0, 1, 0, 0, 1, 0)
  score_residuals <- c(-2, 1, 1, -1, 2, 1)
  corrupt <- list(
    function(x) { x$q <- NULL; x },
    function(x) { x$q <- x$q[-1L, , drop = FALSE]; x },
    function(x) { x$r[1L, 1L] <- 0; x },
    function(x) { x$q[1L, 1L] <- NA_real_; x },
    function(x) { x$z[1L, 1L] <- Inf; x },
    function(x) { x$sqrt_weights[1L] <- 0; x },
    function(x) { x$qr$rank <- 0L; x },
    function(x) { x$qr$pivot <- c(1L, 1L); x },
    function(x) { x$qr$pivot <- c(NA_integer_, 2L); x },
    function(x) { x$qr$pivot <- c(1L, 3L); x }
  )
  for (mutate in corrupt) {
    expect_error(residualized_score_cpp(mutate(projection), exposure,
                                         score_residuals), "cached factors")
  }
})

test_that("fused kernels reject nonfinite products and preserve RNG state", {
  projection <- prepare_projection(matrix(1, 4, 1), rep(1, 4))
  exposure <- c(-1, 1, -1, 1)
  expect_error(residualized_score_cpp(projection, exposure * 1e200,
                                      rep(1e200, 4)), "contributions must be finite")
  set.seed(821)
  before <- .Random.seed
  residualized_score_cpp(projection, exposure, c(1, 2, 3, 4))
  residualized_score_batch_cpp(projection, cbind(exposure, 2 * exposure),
                                c(1, 2, 3, 4))
  expect_identical(.Random.seed, before)
})
