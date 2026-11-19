test_that("permutation p-values use each requested tail and plus one", {
  permuted <- c(-3, -2, 0, 1, 2, 3)
  expect_equal(permutation_p_value(2, permuted), 5 / 7)
  expect_equal(permutation_p_value(2, permuted, "two.sided"), 5 / 7)
  expect_equal(permutation_p_value(2, permuted, "greater"), 3 / 7)
  expect_equal(permutation_p_value(2, permuted, "less"), 6 / 7)
  expect_identical(permutation_extremes(2, permuted, "two.sided"),
                   c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))
  expect_identical(permutation_extremes(2, permuted, "greater"),
                   c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_identical(permutation_extremes(2, permuted, "less"),
                   c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE))
})

test_that("inclusive ties and one-permutation budgets are handled", {
  for (alternative in c("two.sided", "greater", "less")) {
    expect_equal(permutation_p_value(2, rep(2, 5), alternative), 1)
    expect_equal(permutation_p_value(0, 0, alternative), 1)
  }
  expect_equal(permutation_p_value(2, c(2, -2, 1)), 3 / 4)
  expect_equal(permutation_p_value(2, c(2, -2, 1), "greater"), 1 / 2)
  expect_equal(permutation_p_value(2, c(2, -2, 1), "less"), 1)
  expect_equal(permutation_p_value(2, 1, "greater"), 1 / 2)
  expect_equal(permutation_p_value(2, 1, "less"), 1)
})

test_that("permutation p-values attain the finite-budget minimum", {
  permuted <- c(-2, -1, 0, 1, 2)
  minimum <- 1 / (length(permuted) + 1)
  expect_equal(permutation_p_value(3, permuted), minimum)
  expect_equal(permutation_p_value(3, permuted, "greater"), minimum)
  expect_equal(permutation_p_value(-3, permuted, "less"), minimum)
  expect_equal(permutation_p_value(0, permuted), 1)
})

test_that("sign reversal preserves two-sided p-values and swaps one-sided tails", {
  observed <- 1.5
  permuted <- c(-3, -1.5, -0.2, 0.7, 1.5, 2)
  expect_equal(permutation_p_value(observed, permuted),
               permutation_p_value(-observed, -permuted))
  expect_equal(permutation_p_value(observed, permuted, "greater"),
               permutation_p_value(-observed, -permuted, "less"))
  expect_equal(permutation_p_value(observed, permuted, "less"),
               permutation_p_value(-observed, -permuted, "greater"))
})

test_that("numeric names and integer inputs do not change the scalar result", {
  observed <- c(observed = 2L)
  permuted <- c(first = -3L, second = 0L, third = 2L)
  expect_identical(permutation_p_value(observed, permuted), 3 / 4)
  expect_null(names(permutation_p_value(observed, permuted)))
})

test_that("invalid observed and permutation statistics are rejected", {
  invalid_observed <- list(
    NULL, numeric(), c(1, 2), NA_real_, NaN, Inf, -Inf, 1i, TRUE,
    "1", list(1), matrix(1), factor("a"), as.Date("2026-01-01")
  )
  for (observed in invalid_observed) {
    expect_error(permutation_p_value(observed, c(0, 1)), "observed must")
  }
  invalid_permuted <- list(
    NULL, numeric(), c(1, NA_real_), c(1L, NA_integer_), c(1, NaN),
    c(1, Inf), c(1, -Inf), c(1i, 2i), c(TRUE, FALSE), c("1", "2"),
    list(1, 2), matrix(1:4, 2), array(1:4, c(2, 1, 2)),
    factor(c("a", "b")), as.Date(c("2026-01-01", "2026-01-02"))
  )
  for (permuted in invalid_permuted) {
    expect_error(permutation_p_value(1, permuted), "permuted must")
  }
  for (alternative in list("invalid", NA_character_, 1, c("greater", "less"))) {
    expect_error(permutation_p_value(1, c(0, 1), alternative))
  }
})

test_that("permutation p-value calculation does not consume random numbers", {
  set.seed(658)
  before <- .Random.seed
  for (alternative in c("two.sided", "greater", "less")) {
    permutation_p_value(2, c(-3, 0, 2), alternative)
  }
  expect_identical(.Random.seed, before)
})
