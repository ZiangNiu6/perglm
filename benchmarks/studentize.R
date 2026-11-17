# Run with an installed development version of perglm:
# Rscript benchmarks/studentize.R
# Timings exclude compilation, input generation, and correctness checks.
# Allocations are R-heap allocations, not a measurement of native peak memory.

library(perglm)

studentize_r <- getFromNamespace("studentize_score", "perglm")
studentize_cpp <- getFromNamespace("studentize_score_cpp", "perglm")
studentize_batch <- getFromNamespace("studentize_score_batch_cpp", "perglm")

median_elapsed <- function(fun, repetitions = 5L) {
  median(replicate(repetitions, {
    gc()
    unname(system.time(fun())["elapsed"])
  }))
}

allocated_bytes <- function(fun) {
  if (!capabilities("profmem")) {
    return(NA_real_)
  }
  profile <- tempfile("perglm-allocations-")
  on.exit(unlink(profile))
  Rprofmem(profile)
  tryCatch(invisible(fun()), finally = Rprofmem(NULL))
  entries <- readLines(profile)
  entries <- entries[grepl("^[0-9]+ ", entries)]
  sum(as.numeric(sub(" .*", "", entries)))
}

benchmark_case <- function(n, columns) {
  contributions <- matrix(rnorm(n * columns), nrow = n)
  vectors <- lapply(seq_len(columns), function(column) contributions[, column])
  r_calls <- function() vapply(vectors, studentize_r, numeric(1))
  cpp_calls <- function() vapply(vectors, studentize_cpp, numeric(1))
  cpp_batch <- function() studentize_batch(contributions)

  reference <- r_calls()
  stopifnot(isTRUE(all.equal(reference, cpp_calls(), tolerance = 1e-12)))
  stopifnot(isTRUE(all.equal(reference, cpp_batch(), tolerance = 1e-12)))

  data.frame(
    n = n,
    columns = columns,
    r_seconds = median_elapsed(r_calls),
    cpp_calls_seconds = median_elapsed(cpp_calls),
    cpp_batch_seconds = median_elapsed(cpp_batch),
    r_heap_bytes = allocated_bytes(r_calls),
    cpp_calls_r_heap_bytes = allocated_bytes(cpp_calls),
    cpp_batch_r_heap_bytes = allocated_bytes(cpp_batch),
    max_absolute_error = max(abs(reference - cpp_batch()))
  )
}

set.seed(901)
results <- rbind(
  benchmark_case(100L, 5000L),
  benchmark_case(1000L, 1000L),
  benchmark_case(10000L, 200L)
)
print(results, row.names = FALSE)

# Measure the complete projection -> contribution -> studentization path.
prepare <- getFromNamespace("prepare_projection", "perglm")
residualize <- getFromNamespace("residualize_exposure", "perglm")
score_cpp <- getFromNamespace("residualized_score_cpp", "perglm")
score_batch <- getFromNamespace("residualized_score_batch_cpp", "perglm")

benchmark_projection <- function(n, p, columns) {
  z <- cbind(1, matrix(rnorm(n * (p - 1L)), nrow = n))
  projection <- prepare(z, exp(rnorm(n, sd = 0.3)))
  score_residuals <- rnorm(n)
  exposures <- matrix(rnorm(n * columns), nrow = n)
  vectors <- lapply(seq_len(columns), function(column) exposures[, column])
  r_calls <- function() vapply(vectors, function(exposure) {
    studentize_r(residualize(projection, exposure) * score_residuals)
  }, numeric(1))
  cpp_calls <- function() vapply(vectors, function(exposure) {
    score_cpp(projection, exposure, score_residuals)
  }, numeric(1))
  cpp_batch <- function() score_batch(projection, exposures, score_residuals)

  reference <- r_calls()
  stopifnot(isTRUE(all.equal(reference, cpp_calls(), tolerance = 1e-11)))
  stopifnot(isTRUE(all.equal(reference, cpp_batch(), tolerance = 1e-11)))

  data.frame(
    n = n,
    p = p,
    columns = columns,
    r_seconds = median_elapsed(r_calls),
    cpp_calls_seconds = median_elapsed(cpp_calls),
    cpp_batch_seconds = median_elapsed(cpp_batch),
    r_heap_bytes = allocated_bytes(r_calls),
    cpp_calls_r_heap_bytes = allocated_bytes(cpp_calls),
    cpp_batch_r_heap_bytes = allocated_bytes(cpp_batch),
    max_absolute_error = max(abs(reference - cpp_batch()))
  )
}

projection_results <- rbind(
  benchmark_projection(1000L, 5L, 200L),
  benchmark_projection(5000L, 8L, 200L),
  benchmark_projection(20000L, 8L, 100L)
)
print(projection_results, row.names = FALSE)
