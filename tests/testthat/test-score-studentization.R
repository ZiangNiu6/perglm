test_that("studentization uses the empirical variance with divisor n", {
  expect_equal(studentize_score(c(1, 2, 3)), 6 / sqrt(2), tolerance = 1e-14)
  expect_equal(studentize_score(c(-2, -1, 1, 2)), 0)
})

test_that("studentization is stable under extreme common scaling", {
  contributions <- c(-2, 1, 3, 5)
  reference <- sum(contributions) /
    sqrt(sum((contributions - mean(contributions))^2))
  for (scale in c(1e-250, 1, 1e250)) {
    expect_equal(studentize_score(scale * contributions), reference,
                 tolerance = 1e-14)
  }
  expect_equal(studentize_score(-contributions), -reference, tolerance = 1e-14)
})

test_that("studentization preserves tiny differences between large or tiny scores", {
  smallest <- .Machine$double.xmin * .Machine$double.eps
  for (contributions in list(c(1e308, 1e308 + 2e292),
                            c(1e-308, 1e-308 + smallest))) {
    scale <- max(abs(contributions))
    # For two scores, the centered norm is abs(r[2] - r[1]) / sqrt(2).
    expected <- sqrt(2) * sum(contributions / scale) /
      (abs(contributions[2] - contributions[1]) / scale)
    expect_equal(studentize_score(contributions), expected, tolerance = 1e-14)
  }
})

test_that("studentization handles anchor differences that overflow", {
  contributions <- c(1e308, -1e308, 1e308)
  expect_equal(studentize_score(contributions), sqrt(3 / 8), tolerance = 1e-14)
})

test_that("studentization rejects invalid and degenerate contributions", {
  for (contributions in list(numeric(), 1, c(1, NA), c(1, Inf),
                            c(1, NaN), c(1i, 2i), matrix(1:4, 2, 2))) {
    expect_error(studentize_score(contributions), "contributions must")
  }
  for (contributions in list(c(0, 0), rep(2, 4), rep(-2, 4))) {
    expect_error(studentize_score(contributions), "positive empirical variance")
  }
})
