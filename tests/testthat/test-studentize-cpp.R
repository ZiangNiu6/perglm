test_that("the compiled studentizer agrees with the R reference", {
  smallest <- .Machine$double.xmin * .Machine$double.eps
  cases <- list(
    c(1L, 2L, 3L),
    c(a = 1, b = 2, c = 3),
    c(-2, -1, 1, 2),
    c(1e308, 1e308 + 2e292),
    c(1e-308, 1e-308 + smallest),
    c(-1.7e308, 1.7e308, 0),
    smallest * c(1, 2, 3),
    c(1e16, 1, -1e16),
    c(1e16, -1e16, 1)
  )
  for (contributions in cases) {
    reference <- studentize_score(contributions)
    compiled <- studentize_score_cpp(contributions)
    # Compensated C++ sums can differ from R by rounding near cancellation.
    expect_lte(abs(compiled - reference), 1e-13 * max(1, abs(reference)))
  }
  expect_equal(studentize_score_cpp(c(1, 2, 3)), 6 / sqrt(2), tolerance = 1e-14)
})

test_that("the compiled studentizer preserves scale and sign properties", {
  contributions <- c(-2, 1, 3, 5)
  reference <- studentize_score(contributions)
  for (scale in c(1e-250, 1, 1e250)) {
    expect_equal(studentize_score_cpp(scale * contributions), reference,
                 tolerance = 1e-14)
  }
  expect_equal(studentize_score_cpp(-contributions), -reference,
               tolerance = 1e-14)
})

test_that("compiled and R studentizers reject the same invalid inputs", {
  invalid <- list(
    numeric(), 1, c(1, NA_real_), c(1L, NA_integer_), c(1, Inf),
    c(1, NaN), c(1i, 2i), c(TRUE, FALSE), c("1", "2"),
    as.raw(1:3), list(1, 2), matrix(1:4, 2, 2), factor(c("a", "b")),
    as.Date(c("2026-01-01", "2026-01-02")),
    as.POSIXct(c("2026-01-01", "2026-01-02"), tz = "UTC"),
    as.difftime(c(1, 2), units = "days")
  )
  for (contributions in invalid) {
    expect_error(studentize_score(contributions), "contributions must")
    expect_error(studentize_score_cpp(contributions), "contributions must")
  }
  for (contributions in list(c(0, 0), rep(2, 4), rep(-2, 4))) {
    expect_error(studentize_score_cpp(contributions), "positive empirical variance")
  }
})

test_that("the batch studentizer agrees with individual calls and R", {
  set.seed(281)
  contributions <- matrix(rnorm(300), nrow = 30)
  colnames(contributions) <- paste0("score", seq_len(ncol(contributions)))
  original <- contributions
  reference <- apply(contributions, 2L, studentize_score)
  individual <- vapply(seq_len(ncol(contributions)), function(column) {
    studentize_score_cpp(contributions[, column])
  }, numeric(1))
  batch <- studentize_score_batch_cpp(contributions)
  expect_equal(batch, reference, tolerance = 1e-13)
  expect_equal(unname(batch), individual, tolerance = 1e-14)
  expect_identical(contributions, original)
  expect_identical(studentize_score_batch_cpp(matrix(numeric(), 2, 0)), numeric())
  expect_equal(studentize_score_batch_cpp(matrix(c(1, 2, 3), 3, 1)),
               studentize_score_cpp(c(1, 2, 3)))

  integer_matrix <- matrix(1:12, 3, 4)
  expect_equal(studentize_score_batch_cpp(integer_matrix),
               apply(integer_matrix, 2L, studentize_score), tolerance = 1e-13)
})

test_that("the batch studentizer handles extreme columns", {
  contributions <- cbind(
    c(1, 2, 3),
    c(1, 2, 3) * 1e250,
    c(1, 2, 3) * 1e-250,
    c(1e308, -1e308, 1e308)
  )
  expect_equal(studentize_score_batch_cpp(contributions),
               apply(contributions, 2L, studentize_score), tolerance = 1e-13)
})

test_that("the batch studentizer rejects invalid or degenerate columns", {
  for (contributions in list(1:4, matrix(1:4, 1, 4), matrix(1:4 + 0i, 2, 2),
                            matrix(TRUE, 2, 2), matrix(c(1L, NA_integer_), 2),
                            matrix(c(1, Inf), 2))) {
    expect_error(studentize_score_batch_cpp(contributions), "contributions must")
  }
  expect_error(studentize_score_batch_cpp(cbind(c(1, 2, 3), c(0, 0, 0))),
               "Column 2:.*positive empirical variance")
})

test_that("compiled studentization does not consume random numbers", {
  set.seed(457)
  before <- .Random.seed
  studentize_score_cpp(c(1, 2, 3))
  studentize_score_batch_cpp(matrix(1:6, 3, 2))
  expect_identical(.Random.seed, before)
})
