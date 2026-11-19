#!/usr/bin/env Rscript
# Reproducible statistical validation of the exported perglm_test() interface.
# Run from the package root, using an installed development version:
# R_LIBS=/path/to/library Rscript validation/validate-perglm.R \
#   --repetitions=1000 --B=499 --workers=2 --output=validation/results
# Optional: --case=gaussian_mi_n100_gamma0 (comma-separated IDs are accepted).
# The experiment records every requested dataset, including errors and warnings.
# Simulations describe rejection probabilities; their intervals are not unit tests.

parse_options <- function(arguments) {
  options <- list(repetitions = 1000L, B = 499L, workers = 2L,
                  output = "validation/results", case = NULL)
  if ("--help" %in% arguments) {
    cat("Usage: Rscript validation/validate-perglm.R [--repetitions=1000]",
        "[--B=499] [--workers=2] [--output=validation/results] [--case=ID,...]\n")
    quit(status = 0L)
  }
  for (argument in arguments) {
    position <- regexpr("=", argument, fixed = TRUE)[1L]
    if (!startsWith(argument, "--") || position < 4L) {
      stop("Arguments must use --name=value syntax.", call. = FALSE)
    }
    name <- substr(argument, 3L, position - 1L)
    value <- substring(argument, position + 1L)
    if (!name %in% names(options) || !nzchar(value)) {
      stop("Unknown or empty option: ", argument, call. = FALSE)
    }
    if (name %in% c("repetitions", "B", "workers")) {
      number <- suppressWarnings(as.numeric(value))
      maximum <- if (name == "repetitions") 1000000 else .Machine$integer.max
      if (length(number) != 1L || !is.finite(number) || number < 1 ||
          number > maximum || number != floor(number)) {
        stop(name, " must be a positive integer at most ", maximum, ".",
             call. = FALSE)
      }
      options[[name]] <- as.integer(number)
    } else {
      options[[name]] <- value
    }
  }
  options
}

make_cases <- function() {
  definitions <- data.frame(
    scenario = c("gaussian_mi", "gaussian_cm", "gaussian_cm_hetero",
                 "poisson_mi", "poisson_cm", "poisson_cm_nb"),
    family = c(rep("gaussian", 3L), rep("poisson", 3L)),
    target = c("MI", "CM", "CM_only", "MI", "CM", "CM"),
    assignment = c("independent", "confounded", "confounded",
                   "independent", "confounded", "confounded"),
    mean_specification = c("omitted_quadratic", "correct", "correct",
                           "omitted_quadratic", "correct", "correct"),
    outcome_distribution = c("Gaussian", "Gaussian", "heteroskedastic_Gaussian",
                             "Poisson", "Poisson", "negative_binomial_size2"),
    stringsAsFactors = FALSE
  )
  nulls <- do.call(rbind, lapply(seq_len(nrow(definitions)), function(i) {
    result <- definitions[rep(i, 2L), , drop = FALSE]
    result$n <- c(100L, 500L)
    result$gamma <- 0
    result$category <- "regular_null"
    result
  }))
  power <- do.call(rbind, lapply(c("gaussian_cm", "poisson_cm"), function(scenario) {
    effects <- if (scenario == "gaussian_cm") c(0.25, 0.5) else c(0.3, 0.6)
    grid <- expand.grid(n = c(100L, 500L), gamma = effects)
    result <- definitions[rep(match(scenario, definitions$scenario), nrow(grid)), ,
                          drop = FALSE]
    result$target <- "alternative"
    result$n <- grid$n
    result$gamma <- grid$gamma
    result$category <- "power"
    result
  }))
  stress <- definitions[c(1L, 4L), , drop = FALSE]
  stress$scenario <- c("gaussian_rare_x", "poisson_rare_x")
  stress$target <- "MI_stress"
  stress$assignment <- "rare_independent"
  stress$mean_specification <- "correct"
  stress$n <- 40L
  stress$gamma <- 0
  stress$category <- "stress"
  cases <- rbind(nulls, power, stress)
  rownames(cases) <- NULL
  cases$case_index <- seq_len(nrow(cases))
  cases$case_id <- paste0(cases$scenario, "_n", cases$n, "_gamma", cases$gamma)
  cases[, c("case_id", "case_index", "category", "scenario", "family", "n",
             "gamma", "target", "assignment", "mean_specification",
             "outcome_distribution")]
}

generate_data <- function(case) {
  n <- case$n
  u <- stats::runif(n, -1, 1)
  v <- stats::runif(n, -1, 1)
  probability <- switch(case$assignment,
    independent = rep(0.3, n),
    confounded = stats::plogis(stats::qlogis(0.3) + 1.2 * u - 0.8 * v),
    rare_independent = rep(0.02, n)
  )
  x <- stats::rbinom(n, 1L, probability)
  quadratic <- if (case$mean_specification == "omitted_quadratic") 0.6 * u^2 else 0
  if (case$family == "gaussian") {
    mean <- 0.4 + 0.6 * u - 0.4 * v + case$gamma * x + quadratic
    sd <- if (case$scenario == "gaussian_cm_hetero") exp(0.5 * u + 0.3 * x) else 1
    y <- mean + sd * stats::rnorm(n)
  } else {
    mean <- exp(0.2 + 0.5 * u - 0.3 * v + case$gamma * x + quadratic)
    y <- if (case$scenario == "poisson_cm_nb") {
      stats::rnbinom(n, mu = mean, size = 2)
    } else {
      stats::rpois(n, mean)
    }
  }
  data.frame(y = y, u = u, v = v, x = x)
}

run_case <- function(case, repetitions, B, base_seed, alpha) {
  family <- if (case$family == "gaussian") stats::gaussian() else stats::poisson()
  rows <- vector("list", repetitions)
  started <- proc.time()[["elapsed"]]
  for (replicate in seq_len(repetitions)) {
    # Each case has a reserved block covering the maximum supported 1e6 draws.
    # Two adjacent seeds keep data generation separate from permutation sampling.
    data_seed <- as.integer(base_seed + (case$case_index - 1) * 2000002 +
                              (replicate - 1) * 2)
    permutation_seed <- data_seed + 1L
    warnings <- character()
    exposed <- NA_integer_
    result <- tryCatch(withCallingHandlers({
      set.seed(data_seed)
      data <- generate_data(case)
      exposed <- sum(data$x)
      set.seed(permutation_seed)
      test <- perglm::perglm_test(
        y ~ u + v, data = data, exposure = data$x, family = family,
        B = B, alternative = "two.sided", batch_size = 64L,
        keep_statistics = FALSE, control = list(epsilon = 1e-10, maxit = 100L)
      )
      if (!inherits(test, "htest") || length(test$p.value) != 1L ||
          !is.finite(test$p.value) || test$p.value < 1 / (B + 1) ||
          test$p.value > 1 || !is.finite(test$statistic)) {
        stop("The public result violates the finite htest/p-value contract.")
      }
      test
    }, warning = function(warning) {
      warnings <<- c(warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }), error = function(error) error)
    ok <- !inherits(result, "error")
    rows[[replicate]] <- data.frame(
      case_id = case$case_id, replicate = replicate,
      data_seed = data_seed, permutation_seed = permutation_seed,
      exposed = exposed, status = if (ok) "ok" else "error",
      p_value = if (ok) result$p.value else NA_real_,
      statistic = if (ok) unname(result$statistic) else NA_real_,
      exceedances = if (ok) result$exceedances else NA_real_,
      rejected = if (ok) result$p.value <= alpha else FALSE,
      nuisance_score = if (ok) result$null_diagnostics$max_score_correlation else NA_real_,
      fit_iterations = if (ok) result$null_diagnostics$iterations else NA_integer_,
      warning_count = length(warnings),
      warnings = paste(unique(warnings), collapse = " | "),
      error = if (ok) "" else conditionMessage(result),
      stringsAsFactors = FALSE
    )
  }
  result <- do.call(rbind, rows)
  elapsed <- proc.time()[["elapsed"]] - started
  cat(sprintf("Completed %s: %d datasets, %d errors, %.2f s\n",
              case$case_id, repetitions, sum(result$status == "error"), elapsed))
  flush.console()
  result
}

summarize_case <- function(case, rows, alpha) {
  requested <- nrow(rows)
  successful <- sum(rows$status == "ok")
  rejected <- sum(rows$rejected)
  interval <- stats::binom.test(rejected, requested, conf.level = 0.95)$conf.int
  data.frame(
    case, requested = requested, successful = successful,
    errors = requested - successful, warning_datasets = sum(rows$warning_count > 0L),
    rejections = rejected, alpha = alpha,
    unconditional_rejection_rate = rejected / requested,
    rejection_ci_lower = unname(interval[1L]),
    rejection_ci_upper = unname(interval[2L]),
    success_rate = successful / requested,
    success_conditional_rejection_rate = if (successful) rejected / successful else NA_real_,
    mean_successful_p_value = if (successful) mean(rows$p_value[rows$status == "ok"]) else NA_real_,
    stringsAsFactors = FALSE
  )
}

options <- parse_options(commandArgs(trailingOnly = TRUE))
cases <- make_cases()
if (!is.null(options$case)) {
  selected <- strsplit(options$case, ",", fixed = TRUE)[[1L]]
  if (any(!selected %in% cases$case_id)) {
    stop("Unknown case ID(s): ", paste(setdiff(selected, cases$case_id), collapse = ", "),
         "\nValid IDs:\n", paste(cases$case_id, collapse = "\n"), call. = FALSE)
  }
  cases <- cases[cases$case_id %in% selected, , drop = FALSE]
}
if (!requireNamespace("perglm", quietly = TRUE) ||
    !"perglm_test" %in% getNamespaceExports("perglm")) {
  stop("Install the development package exporting perglm_test() before running.",
       call. = FALSE)
}
effective_workers <- if (.Platform$OS.type == "windows") 1L else
  min(options$workers, nrow(cases))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
base_seed <- 20261118L
alpha <- 0.05
dir.create(options$output, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(options$output)) {
  stop("Cannot create output directory: ", options$output, call. = FALSE)
}
script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[1L]), mustWork = TRUE)
project_root <- dirname(dirname(script_path))
source_paths <- c(list.files(file.path(project_root, "R"), pattern = "\\.R$", full.names = TRUE),
                  list.files(file.path(project_root, "src"),
                             pattern = "\\.(cpp|h)$|^Makevars(\\.win)?$", full.names = TRUE),
                  script_path, file.path(project_root, "DESCRIPTION"),
                  file.path(project_root, "NAMESPACE"))
source_hashes <- tools::md5sum(source_paths)
names(source_hashes) <- substring(names(source_hashes), nchar(project_root) + 2L)
config <- list(
  requested_options = options, effective_workers = effective_workers,
  alpha = alpha, base_seed = base_seed, rng_kind = RNGkind(),
  seed_rule = "data = base_seed + (case_index-1)*2000002 + (replicate-1)*2; permutation=data+1",
  interval = "95% pointwise exact binomial (Clopper-Pearson); no pass/fail rule",
  failure_rule = "Every requested dataset retained; an error gives no decision and no rejection",
  R_version = R.version.string, platform = R.version$platform,
  perglm_version = as.character(utils::packageVersion("perglm")),
  installed_package_path = normalizePath(find.package("perglm")),
  source_md5 = source_hashes,
  started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
saveRDS(config, file.path(options$output, "config.rds"))
utils::write.csv(cases, file.path(options$output, "scenarios.csv"), row.names = FALSE)
capture.output(utils::sessionInfo(), file = file.path(options$output, "session-info.txt"))
cat(sprintf("Validating %d cases x %d datasets; B=%d; workers=%d\n",
            nrow(cases), options$repetitions, options$B, effective_workers))
flush.console()
job <- function(index) {
  run_case(cases[index, , drop = FALSE], options$repetitions, options$B, base_seed, alpha)
}
started <- proc.time()[["elapsed"]]
results <- if (effective_workers > 1L) {
  parallel::mclapply(seq_len(nrow(cases)), job, mc.cores = effective_workers,
                     mc.preschedule = FALSE, mc.set.seed = FALSE)
} else {
  lapply(seq_len(nrow(cases)), job)
}
if (any(vapply(results, inherits, logical(1), what = "try-error"))) {
  stop("A simulation worker failed unexpectedly; results were not summarized.", call. = FALSE)
}
replicates <- do.call(rbind, results)
rownames(replicates) <- NULL
summary <- do.call(rbind, lapply(seq_len(nrow(cases)), function(index) {
  summarize_case(cases[index, , drop = FALSE], results[[index]], alpha)
}))
rownames(summary) <- NULL
connection <- gzfile(file.path(options$output, "replicates.csv.gz"), open = "wt")
tryCatch(utils::write.csv(replicates, connection, row.names = FALSE),
         finally = close(connection))
utils::write.csv(summary, file.path(options$output, "summary.csv"), row.names = FALSE)
config$completed_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
config$elapsed_seconds <- proc.time()[["elapsed"]] - started
saveRDS(config, file.path(options$output, "config.rds"))
cat("\nPointwise intervals describe simulation uncertainty, not proof of calibration.\n")
print(summary[, c("case_id", "successful", "errors", "rejections",
                  "unconditional_rejection_rate", "rejection_ci_lower",
                  "rejection_ci_upper", "success_conditional_rejection_rate")],
      row.names = FALSE, digits = 4)
cat(sprintf("\nFinished in %.2f s. Results: %s\n", config$elapsed_seconds,
            normalizePath(options$output)))
