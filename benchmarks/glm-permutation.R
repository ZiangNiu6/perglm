# Run against an installed development version of perglm:
# R_LIBS=/path/to/library Rscript benchmarks/glm-permutation.R
# Compilation, data generation, null fitting, preparation and correctness checks
# are outside the timed region. Each timed method includes the same seeded
# sample.int() permutation draws, scoring, and p-value calculation.
# Rprofmem measures R-heap allocations, not native allocations or peak memory.

library(perglm)

fit_null <- getFromNamespace("fit_glm_null", "perglm")
residualize_r <- getFromNamespace("residualize_exposure", "perglm")
studentize_r <- getFromNamespace("studentize_score", "perglm")
calibrate <- getFromNamespace("calibrate_glm_permutation", "perglm")

n <- 5000L
p <- 8L
B <- 200L
batch_size <- 64L
repetitions <- 3L
data_seed <- 20261118L
permutation_seed <- 9417L

set.seed(data_seed)
covariates <- matrix(rnorm(n * (p - 1L)), nrow = n)
colnames(covariates) <- paste0("z", seq_len(p - 1L))
mu <- exp(0.3 + 0.2 * covariates[, 1L] - 0.15 * covariates[, 2L])
data <- data.frame(y = rpois(n, mu), covariates)
exposure <- as.numeric(rbinom(n, 1L, 0.1))
null <- fit_null(y ~ ., data, family = stats::poisson(),
                  control = stats::glm.control(epsilon = 1e-10))
stopifnot(ncol(null$projection$z) == p)

# The reference reuses the prepared weighted QR and fixed response residuals.
# It residualizes raw shuffled labels and never refits the outcome model.
r_baseline <- function() {
  set.seed(permutation_seed)
  observed <- studentize_r(
    residualize_r(null$projection, exposure) * null$score_residuals
  )
  statistics <- numeric(B)
  for (draw in seq_len(B)) {
    shuffled <- exposure[sample.int(n, n, replace = FALSE)]
    statistics[draw] <- studentize_r(
      residualize_r(null$projection, shuffled) * null$score_residuals
    )
  }
  exceedances <- sum(abs(statistics) >= abs(observed))
  list(statistic = observed, p.value = (1 + exceedances) / (B + 1),
       exceedances = exceedances, permuted_statistics = statistics,
       rng_state = .Random.seed)
}

native_calibration <- function(keep_statistics) {
  set.seed(permutation_seed)
  result <- calibrate(null, exposure, B = B, alternative = "two.sided",
                       batch_size = batch_size,
                       keep_statistics = keep_statistics)
  result$rng_state <- .Random.seed
  result
}

discard_statistics <- function() native_calibration(FALSE)
retain_statistics <- function() native_calibration(TRUE)

assert_close <- function(actual, expected, tolerance = 1e-11) {
  stopifnot(length(actual) == length(expected),
            all(is.finite(actual)), all(is.finite(expected)),
            all(abs(actual - expected) <= tolerance * pmax(1, abs(expected))))
}

# These initial calls also warm up all paths before timing.
reference <- r_baseline()
discarded <- discard_statistics()
retained <- retain_statistics()
for (result in list(discarded, retained)) {
  assert_close(result$statistic, reference$statistic)
  stopifnot(identical(result$p.value, reference$p.value),
            result$exceedances == reference$exceedances,
            identical(result$rng_state, reference$rng_state))
}
stopifnot(is.null(discarded$permuted_statistics))
assert_close(retained$permuted_statistics, reference$permuted_statistics)

median_elapsed <- function(fun) {
  median(replicate(repetitions, {
    gc()
    unname(system.time(fun())["elapsed"])
  }))
}

allocated_bytes <- function(fun) {
  if (!capabilities("profmem")) {
    return(NA_real_)
  }
  profile <- tempfile("perglm-permutation-allocations-")
  on.exit(unlink(profile))
  Rprofmem(profile)
  tryCatch(invisible(fun()), finally = Rprofmem(NULL))
  entries <- readLines(profile)
  entries <- entries[grepl("^[0-9]+ ", entries)]
  sum(as.numeric(sub(" .*", "", entries)))
}

methods <- list(r_baseline = r_baseline,
                 native_discard_statistics = discard_statistics,
                 native_retain_statistics = retain_statistics)
seconds <- vapply(methods, median_elapsed, numeric(1))
allocations <- vapply(methods, allocated_bytes, numeric(1))
results <- data.frame(
  method = names(methods), n = n, p = p, B = B, batch_size = batch_size,
  repetitions = repetitions, median_seconds = unname(seconds),
  speedup_over_r = unname(seconds[1L] / seconds),
  r_heap_bytes = unname(allocations),
  retains_statistics = c(TRUE, FALSE, TRUE),
  p_value = reference$p.value,
  max_absolute_statistic_error = c(
    0, NA_real_,
    max(abs(retained$permuted_statistics - reference$permuted_statistics))
  )
)

cat(R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("perglm:", as.character(utils::packageVersion("perglm")), "\n")
cat("Data seed:", data_seed, "; permutation seed:", permutation_seed, "\n")
print(results, row.names = FALSE, digits = 6)
