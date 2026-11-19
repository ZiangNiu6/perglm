#!/usr/bin/env Rscript
# Reproduce finite-fit and numerical-degeneracy guard checks.
# Run from the package root against the installed development version:
# R_LIBS=/path/to/library Rscript validation/edge-probes.R
# All public tests use seed 17 and B = 9; maxit is recorded for every case.
# `passed` compares the observed behavior with the stated mathematical contract.
# A FALSE entry records a regression; package unit tests also cover the guards.

library(perglm)

capture_probe <- function(expression) {
  warnings <- character()
  result <- tryCatch(withCallingHandlers(expression, warning = function(warning) {
    warnings <<- c(warnings, conditionMessage(warning))
    invokeRestart("muffleWarning")
  }), error = function(error) error)
  list(result = result, warnings = unique(warnings))
}

edge_cases <- function() {
  small <- data.frame(y = rep(0, 20), z = seq(-1, 1, length.out = 20),
                      x = sin(seq_len(20)))
  add_case <- function(name, data, family, expected, notes,
                       formula = y ~ z, maxit = 25L,
                       family_name = family$family) {
    list(case = name, data = data, formula = formula, family = family,
         family_name = family_name, expected = expected, notes = notes,
         maxit = maxit)
  }
  intercept_error <- "error: no finite null-model MLE"
  variance_error <- "error: zero empirical score variance"
  cases <- list(
    add_case("poisson_all_zero", small, stats::poisson(), intercept_error,
      "Exact all-zero count guard; response has no finite log-link intercept MLE."),
    add_case("binomial_all_zero", small, stats::binomial(), intercept_error,
      "Exact Bernoulli endpoint guard; both response classes are required."),
    add_case("fixed_NB_all_zero", small, MASS::negative.binomial(1.7),
      intercept_error, "Fixed-size NB with an all-zero response; theta = 1.7.",
      family_name = "negative.binomial")
  )
  all_one <- small
  all_one$y[] <- 1
  cases[[length(cases) + 1L]] <- add_case("binomial_all_one", all_one,
    stats::binomial(), intercept_error, "Exact Bernoulli endpoint guard.")
  constant <- small
  constant$y[] <- 3
  constant$off <- rep(0.25, nrow(constant))
  cases[[length(cases) + 1L]] <- add_case("Gaussian_constant_no_offset", constant,
    stats::gaussian(), variance_error,
    "Constant response and an intercept imply exactly zero outcome scores.")
  cases[[length(cases) + 1L]] <- add_case("Gaussian_constant_constant_offset",
    constant, stats::gaussian(), variance_error,
    "A constant known offset is absorbed by the intercept.",
    formula = y ~ z + offset(off))
  constant$y[] <- 1
  cases[[length(cases) + 1L]] <- add_case("poisson_constant_positive", constant,
    stats::poisson(), variance_error,
    "Positive constant counts admit an exact common-mean fit.")
  constant$y[] <- 2
  cases[[length(cases) + 1L]] <- add_case("fixed_NB_constant_positive", constant,
    MASS::negative.binomial(1.7), variance_error,
    "Fixed-size NB admits an exact common-mean fit; theta = 1.7.",
    family_name = "negative.binomial")

  i <- seq_len(40)
  z <- rep(c(-1, 1), each = 20)
  exposure <- sin(i * 0.71) + cos(i * 0.17)
  complete <- data.frame(y = as.numeric(z > 0), z = z, x = exposure)
  quasi <- data.frame(y = c(rep(c(0, 1), 10), rep(1, 20)),
    group = factor(rep(c("mixed", "pure"), each = 20)), x = exposure)
  for (link in c("logit", "probit", "cauchit", "cloglog")) {
    cases[[length(cases) + 1L]] <- add_case(paste0("complete_", link), complete,
      stats::binomial(link), intercept_error,
      paste("Complete separation: direction (intercept = 0, z = 1) has",
            "strictly positive signed margins. The separation guard must reject",
            "even if glm() reports convergence at maxit = 200."),
      maxit = 200L)
    cases[[length(cases) + 1L]] <- add_case(paste0("quasi_factor_", link), quasi,
      stats::binomial(link), intercept_error,
      paste("Quasi separation: pure-cell indicator direction (0, 1) has signed",
            "margin zero on the mixed cell and one on all 20 pure-one rows."),
      formula = y ~ group, maxit = 200L)
  }
  cases[[length(cases) + 1L]] <- add_case("complete_cauchit_default_control",
    complete, stats::binomial("cauchit"), intercept_error,
    "Default maxit = 25 rejects nonconvergence; the separation guard also rejects converged fits at larger maxit.")

  linear <- data.frame(z = seq(-1, 1, length.out = 40), x = exposure)
  linear$y <- 2 + 3 * linear$z
  cases[[length(cases) + 1L]] <- add_case("Gaussian_exact_nonconstant_linear",
    linear, stats::gaussian(), variance_error,
    paste("Exact y = 2 + 3*z has zero mathematical score variance;",
          "floating QR residuals are checked against the prediction-roundoff envelope."))
  grouped_linear <- data.frame(z = z, y = 3 + 2 * z, x = exposure)
  cases[[length(cases) + 1L]] <- add_case("Gaussian_exact_two_groups",
    grouped_linear, stats::gaussian(), variance_error,
    "Exact y = 3 + 2*z; the residual-resolution guard must reject projection roundoff.")

  offset_only <- data.frame(y = rep(0:1, 20), x = exposure)
  offset_only$off <- ifelse(offset_only$y == 1, 2, -2)
  cases[[length(cases) + 1L]] <- add_case("binomial_offset_only_classification",
    offset_only, stats::binomial(), "finite p-value: valid interior fit",
    paste("Known offsets classify all rows, but the sole free intercept has a",
          "finite MLE. Separation witnesses must exclude the fixed offset."),
    formula = y ~ 1 + offset(off))
  cases
}

normalized_score_diagnostic <- function(fit, family) {
  mu <- fit$fitted.values
  epsilon <- (fit$y - mu) * family$mu.eta(fit$linear.predictors) /
    family$variance(mu)
  if (any(!is.finite(epsilon))) return(NA_real_)
  epsilon_scale <- max(abs(epsilon))
  if (epsilon_scale == 0) return(0)
  scaled_epsilon <- epsilon / epsilon_scale
  z <- fit$x
  scales <- apply(abs(z), 2L, max)
  if (any(scales == 0)) return(NA_real_)
  scaled_z <- sweep(z, 2L, scales, "/")
  max(abs(drop(crossprod(scaled_z, scaled_epsilon))) /
        (sqrt(colSums(scaled_z^2)) * sqrt(sum(scaled_epsilon^2))))
}

run_edge_probe <- function(case) {
  diagnostic <- capture_probe(stats::glm(case$formula, data = case$data,
    family = case$family, control = stats::glm.control(maxit = case$maxit),
    x = TRUE, y = TRUE, model = TRUE, na.action = stats::na.fail))
  fit <- diagnostic$result
  fit_ok <- !inherits(fit, "error")
  set.seed(17)
  before <- .Random.seed
  public <- capture_probe(perglm::perglm_test(case$formula, data = case$data,
    exposure = case$data$x, family = case$family, B = 9L,
    control = list(maxit = case$maxit)))
  result <- public$result
  public_ok <- !inherits(result, "error")
  finite_result <- public_ok && is.finite(result$p.value) &&
    is.finite(result$statistic)
  binomial <- fit_ok && case$family_name == "binomial"
  margins <- if (binomial) {
    (2 * fit$y - 1) * drop(fit$x %*% fit$coefficients)
  } else NA_real_
  signed_eta <- if (binomial) {
    (2 * fit$y - 1) * fit$linear.predictors
  } else NA_real_
  warning_messages <- c(
    if (length(diagnostic$warnings)) paste0("glm: ", diagnostic$warnings),
    if (length(public$warnings)) paste0("public: ", public$warnings)
  )
  data.frame(
    case = case$case, family = case$family_name, link = case$family$link,
    expected = case$expected,
    observed = if (finite_result) "returned_p_value" else {
      if (public_ok) "nonfinite_result" else "error"
    },
    passed = if (startsWith(case$expected, "error:")) !public_ok else finite_result,
    p_value = if (public_ok) result$p.value else NA_real_,
    statistic = if (public_ok) unname(result$statistic) else NA_real_,
    converged = if (fit_ok) fit$converged else NA,
    boundary = if (fit_ok) fit$boundary else NA,
    iterations = if (fit_ok) fit$iter else NA_integer_,
    scorediag = if (public_ok) result$null_diagnostics$max_score_correlation else {
      if (fit_ok) normalized_score_diagnostic(fit, case$family) else NA_real_
    },
    min_fitted_mean = if (fit_ok) min(fit$fitted.values) else NA_real_,
    max_fitted_mean = if (fit_ok) max(fit$fitted.values) else NA_real_,
    max_absolute_response_residual = if (fit_ok) {
      max(abs(fit$y - fit$fitted.values))
    } else NA_real_,
    min_signed_design_margin = if (binomial) min(margins) else NA_real_,
    all_signed_eta_positive = if (binomial) all(signed_eta > 0) else NA,
    seed = 17L, B = 9L, maxit = case$maxit,
    rng_unchanged = identical(before, .Random.seed),
    error = if (public_ok) "" else conditionMessage(result),
    warnings = paste(warning_messages, collapse = " | "),
    notes = case$notes, stringsAsFactors = FALSE
  )
}

output_directory <- file.path("validation", "results")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
results <- do.call(rbind, lapply(edge_cases(), run_edge_probe))
utils::write.csv(results, file.path(output_directory, "edge-probes.csv"),
                 row.names = FALSE, na = "")
metadata <- c(
  paste("Installed package:", find.package("perglm")),
  paste("Package version:", as.character(utils::packageVersion("perglm"))),
  paste("RNG kinds:", paste(RNGkind(), collapse = ", ")),
  "Public permutation seed: 17; B: 9; maxit recorded in CSV.",
  "Diagnostic GLM uses the same family, formula, observations, and controls.",
  "scorediag is the public diagnostic when available; otherwise computed independently.",
  capture.output(sessionInfo())
)
writeLines(metadata, file.path(output_directory, "edge-probes-session-info.txt"))
print(results[, c("case", "expected", "observed", "passed", "p_value")],
      row.names = FALSE)
cat("\nRecorded", nrow(results), "cases;", sum(results$passed),
    "match the expected contract and", sum(!results$passed),
    "fail the expected guard contract.\n")
