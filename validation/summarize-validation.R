#!/usr/bin/env Rscript
# Analyze saved validation results without rerunning or discarding any dataset.
# Rscript validation/summarize-validation.R --results=validation/results \
#   --exact=validation/exact-results
# Uses base/recommended R only. Figures are exported separately and embedded
# into report.html so the report can be opened without a web connection.

arguments <- list(results = "validation/results", exact = "validation/exact-results")
for (argument in commandArgs(trailingOnly = TRUE)) {
  match <- regexec("^--(results|exact)=(.+)$", argument)
  parts <- regmatches(argument, match)[[1L]]
  if (length(parts) != 3L) {
    stop("Use --results=PATH or --exact=PATH.", call. = FALSE)
  }
  arguments[[parts[2L]]] <- parts[3L]
}
results_directory <- normalizePath(arguments$results, mustWork = TRUE)
script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[1L]), mustWork = TRUE)
project_root <- dirname(dirname(script_path))
config <- readRDS(file.path(results_directory, "config.rds"))
scenarios <- utils::read.csv(file.path(results_directory, "scenarios.csv"),
                            stringsAsFactors = FALSE)
summary <- utils::read.csv(file.path(results_directory, "summary.csv"),
                          stringsAsFactors = FALSE)
replicates <- utils::read.csv(gzfile(file.path(results_directory, "replicates.csv.gz")),
                             stringsAsFactors = FALSE)
if (is.null(config$completed_utc)) {
  stop("The simulation is not marked complete in config.rds.", call. = FALSE)
}

source_paths <- file.path(project_root, names(config$source_md5))
current_source_paths <- c(
  list.files(file.path(project_root, "R"), pattern = "\\.R$", full.names = TRUE),
  list.files(file.path(project_root, "src"),
             pattern = "\\.(cpp|h)$|^Makevars(\\.win)?$", full.names = TRUE),
  file.path(project_root, "validation", "validate-perglm.R"),
  file.path(project_root, "DESCRIPTION"), file.path(project_root, "NAMESPACE")
)
if (!setequal(source_paths, current_source_paths)) {
  stop("The set of R/native source files differs from the simulation snapshot.",
       call. = FALSE)
}
if (any(!file.exists(source_paths))) {
  stop("A source file recorded by the simulation is missing.", call. = FALSE)
}
current_hashes <- unname(tools::md5sum(source_paths))
changed_sources <- names(config$source_md5)[current_hashes != unname(config$source_md5)]
if (length(changed_sources)) {
  stop("Source files changed since the simulation began: ",
       paste(changed_sources, collapse = ", "), call. = FALSE)
}
if (nrow(scenarios) != 22L || anyDuplicated(scenarios$case_id) ||
    !setequal(summary$case_id, scenarios$case_id) ||
    !setequal(replicates$case_id, scenarios$case_id)) {
  stop("Expected the complete set of 22 unique validation cases.", call. = FALSE)
}
summary <- summary[match(scenarios$case_id, summary$case_id), , drop = FALSE]
rownames(summary) <- NULL

equal_numeric <- function(actual, expected, description, tolerance = 1e-12) {
  if (length(actual) != length(expected) ||
      any(is.na(actual) != is.na(expected)) ||
      any(abs(actual - expected) > tolerance * pmax(1, abs(expected)), na.rm = TRUE)) {
    stop("Inconsistent saved results: ", description, call. = FALSE)
  }
}
for (index in seq_len(nrow(scenarios))) {
  scenario <- scenarios[index, , drop = FALSE]
  rows <- replicates[replicates$case_id == scenario$case_id, , drop = FALSE]
  rows <- rows[order(rows$replicate), , drop = FALSE]
  if (!identical(rows$replicate, seq_len(config$requested_options$repetitions)) ||
      any(!rows$status %in% c("ok", "error")) || anyNA(rows$rejected)) {
    stop("Missing, repeated, or invalid replicate rows for ", scenario$case_id,
         call. = FALSE)
  }
  ok <- rows$status == "ok"
  if (any(!is.finite(rows$p_value[ok])) || any(!is.finite(rows$statistic[ok])) ||
      any(!is.na(rows$p_value[!ok])) || any(rows$rejected[!ok]) ||
      any(rows$rejected[ok] != (rows$p_value[ok] <= config$alpha))) {
    stop("Invalid success/error/rejection records for ", scenario$case_id,
         call. = FALSE)
  }
  if (any(rows$exceedances[ok] < 0) ||
      any(rows$exceedances[ok] > config$requested_options$B) ||
      any(rows$exceedances[ok] != floor(rows$exceedances[ok]))) {
    stop("Invalid permutation counts for ", scenario$case_id, call. = FALSE)
  }
  equal_numeric(rows$p_value[ok],
    (1 + rows$exceedances[ok]) / (config$requested_options$B + 1), "plus-one p-values")
  expected_seeds <- config$base_seed + (scenario$case_index - 1) * 2000002 +
    (rows$replicate - 1) * 2
  equal_numeric(rows$data_seed, expected_seeds, "data seeds", tolerance = 0)
  equal_numeric(rows$permutation_seed, expected_seeds + 1, "permutation seeds", tolerance = 0)
  counts <- c(nrow(rows), sum(ok), sum(!ok), sum(rows$warning_count > 0L),
              sum(rows$rejected))
  equal_numeric(unlist(summary[index, c("requested", "successful", "errors",
                                        "warning_datasets", "rejections")]),
                counts, "summary counts", tolerance = 0)
  interval <- stats::binom.test(counts[5L], counts[1L], conf.level = 0.95)$conf.int
  equal_numeric(summary$unconditional_rejection_rate[index], counts[5L] / counts[1L],
                "unconditional rejection rate")
  equal_numeric(summary$success_rate[index], counts[2L] / counts[1L], "success rate")
  equal_numeric(summary$success_conditional_rejection_rate[index],
                if (counts[2L]) counts[5L] / counts[2L] else NA_real_,
                "success-conditional rejection rate")
  equal_numeric(unlist(summary[index, c("rejection_ci_lower", "rejection_ci_upper")]),
                unname(interval), "exact binomial intervals")
  if (!identical(as.character(unlist(summary[index, names(scenarios)])),
                 as.character(unlist(scenario)))) {
    stop("Scenario metadata differ between files for ", scenario$case_id,
         call. = FALSE)
  }
}

nulls <- summary[summary$category == "regular_null", , drop = FALSE]
power <- summary[summary$category == "power", , drop = FALSE]
stress <- summary[summary$category == "stress", , drop = FALSE]
if (nrow(nulls) != 12L || nrow(power) != 8L || nrow(stress) != 2L) {
  stop("Expected 12 regular null, eight power, and two stress cells.", call. = FALSE)
}
nulls$one_sided_excess_p_value <- stats::pbinom(
  nulls$rejections - 1, nulls$requested, config$alpha, lower.tail = FALSE
)
nulls$bonferroni_excess_p_value <- pmin(1, nrow(nulls) * nulls$one_sided_excess_p_value)
nulls$excess_flag <- nulls$bonferroni_excess_p_value < 0.05
utils::write.csv(nulls, file.path(results_directory, "null-calibration-checks.csv"),
                 row.names = FALSE)

scenario_labels <- c(
  gaussian_mi = "MI; omitted quadratic mean",
  gaussian_cm = "CM; Gaussian errors",
  gaussian_cm_hetero = "CM; heteroskedastic errors",
  poisson_mi = "MI; omitted quadratic log mean",
  poisson_cm = "CM; Poisson outcomes",
  poisson_cm_nb = "CM; NB outcomes, Poisson fit",
  gaussian_rare_x = "Rare binary X; Gaussian",
  poisson_rare_x = "Rare binary X; Poisson"
)
series_colors <- c("100" = "#2166AC", "500" = "#B35806")
open_png <- function(filename, width = 1800L, height = 1000L) {
  options <- list(filename = filename, width = width, height = height,
                  res = 180, bg = "white")
  if (capabilities("cairo")) options$type <- "cairo"
  do.call(grDevices::png, options)
}
null_plot <- file.path(results_directory, "null-calibration.png")
open_png(null_plot)
graphics::par(mfrow = c(1, 2), mar = c(4.8, 11.3, 3.5, 1.2),
              oma = c(1.2, 0, 1, 0), cex = 0.85, family = "sans")
for (working_family in c("gaussian", "poisson")) {
  panel <- nulls[nulls$family == working_family, , drop = FALSE]
  positions <- rev(seq_len(nrow(panel)))
  limit <- max(0.085, panel$rejection_ci_upper + 0.005)
  graphics::plot(panel$unconditional_rejection_rate, positions, type = "n",
    xlim = c(0, limit), ylim = c(0.5, nrow(panel) + 0.5), yaxt = "n", ylab = "",
    xlab = "Rejection probability", main = paste0(tools::toTitleCase(working_family),
                                                  " working model"), bty = "n")
  graphics::abline(v = config$alpha, lty = 2, col = "#555555", lwd = 1.2)
  colors <- unname(series_colors[as.character(panel$n)])
  graphics::segments(panel$rejection_ci_lower, positions,
                     panel$rejection_ci_upper, positions, col = colors, lwd = 1.8)
  graphics::segments(panel$rejection_ci_lower, positions - 0.07,
                     panel$rejection_ci_lower, positions + 0.07, col = colors)
  graphics::segments(panel$rejection_ci_upper, positions - 0.07,
                     panel$rejection_ci_upper, positions + 0.07, col = colors)
  graphics::points(panel$unconditional_rejection_rate, positions,
                   pch = ifelse(panel$n == 100L, 16L, 17L), col = colors, cex = 1.1)
  graphics::axis(2, at = positions,
    labels = paste0(unname(scenario_labels[panel$scenario]), "\n(n = ", panel$n, ")"),
    las = 1, tick = FALSE, cex.axis = 0.8)
}
graphics::mtext("95% pointwise exact binomial intervals; dashed line: nominal 0.05",
                side = 1, outer = TRUE, line = 0, cex = 0.82)
grDevices::dev.off()

power_plot <- file.path(results_directory, "power.png")
open_png(power_plot)
graphics::par(mfrow = c(1, 2), mar = c(4.8, 4.8, 3.5, 1.2),
              oma = c(1.2, 0, 1, 0), cex = 0.95, family = "sans")
power_with_baseline <- rbind(
  summary[summary$scenario %in% c("gaussian_cm", "poisson_cm") &
            summary$category == "regular_null", ], power
)
for (working_family in c("gaussian", "poisson")) {
  panel <- power_with_baseline[power_with_baseline$family == working_family, ]
  graphics::plot(NA_real_, NA_real_, type = "n", xlim = c(0, max(panel$gamma) * 1.04),
    ylim = c(0, 1), xlab = "Exposure coefficient on the link scale",
    ylab = "Rejection probability", main = paste0(tools::toTitleCase(working_family),
                                                  " working model"), bty = "l")
  graphics::abline(h = config$alpha, col = "#777777", lty = 3)
  for (n in c(100L, 500L)) {
    series <- panel[panel$n == n, ]
    series <- series[order(series$gamma), ]
    color <- unname(series_colors[as.character(n)])
    graphics::segments(series$gamma, series$rejection_ci_lower,
                       series$gamma, series$rejection_ci_upper, col = color, lwd = 1.5)
    graphics::arrows(series$gamma, series$rejection_ci_lower,
                     series$gamma, series$rejection_ci_upper,
                     angle = 90, code = 3, length = 0.035, col = color)
    graphics::lines(series$gamma, series$unconditional_rejection_rate,
                    type = "b", pch = if (n == 100L) 16L else 17L,
                    lty = if (n == 100L) 1L else 2L, col = color, lwd = 1.8)
  }
  graphics::legend("topleft", legend = c("n = 100", "n = 500"),
                    col = unname(series_colors), pch = c(16L, 17L),
                    lty = c(1L, 2L), bty = "n", cex = 0.85)
}
graphics::mtext("Confounded exposure; correctly specified mean and variance; 95% pointwise intervals",
                side = 1, outer = TRUE, line = 0, cex = 0.82)
grDevices::dev.off()

escape_html <- function(values) {
  values <- as.character(values)
  values[is.na(values)] <- "NA"
  for (pair in list(c("&", "&amp;"), c("<", "&lt;"), c(">", "&gt;"),
                     c('"', "&quot;"), c("'", "&#39;"))) {
    values <- gsub(pair[1L], pair[2L], values, fixed = TRUE)
  }
  values
}
html_table <- function(data) {
  header <- paste0("<th>", escape_html(names(data)), "</th>", collapse = "")
  rows <- vapply(seq_len(nrow(data)), function(index) {
    paste0("<tr>", paste0("<td>", escape_html(unlist(data[index, ], use.names = FALSE)),
                          "</td>", collapse = ""), "</tr>")
  }, character(1))
  paste0('<div class="table-wrap"><table><thead><tr>', header,
         "</tr></thead><tbody>", paste0(rows, collapse = "\n"),
         "</tbody></table></div>")
}
percentage <- function(x) ifelse(is.na(x), "NA", sprintf("%.1f%%", 100 * x))
interval_text <- function(lower, upper) paste0("[", percentage(lower), ", ",
                                              percentage(upper), "]")
format_p <- function(x) ifelse(x < 0.0001, "<0.0001", sprintf("%.4f", x))
png_base64 <- function(path) {
  bytes <- as.integer(readBin(path, what = "raw", n = file.info(path)$size))
  padding <- (3L - length(bytes) %% 3L) %% 3L
  groups <- matrix(c(bytes, rep(0L, padding)), ncol = 3L, byrow = TRUE)
  indexes <- cbind(bitwShiftR(groups[, 1L], 2L),
    bitwOr(bitwShiftL(bitwAnd(groups[, 1L], 3L), 4L), bitwShiftR(groups[, 2L], 4L)),
    bitwOr(bitwShiftL(bitwAnd(groups[, 2L], 15L), 2L), bitwShiftR(groups[, 3L], 6L)),
    bitwAnd(groups[, 3L], 63L))
  alphabet <- strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/", "")[[1L]]
  encoded <- alphabet[as.vector(t(indexes)) + 1L]
  if (padding) encoded[seq.int(length(encoded) - padding + 1L, length(encoded))] <- "="
  paste0(encoded, collapse = "")
}
embedded_image <- function(path, alternative) {
  paste0('<img src="data:image/png;base64,', png_base64(path), '" alt="',
         escape_html(alternative), '">')
}
test_summary_path <- file.path(results_directory, "test-summary.csv")
test_summary <- if (file.exists(test_summary_path)) utils::read.csv(test_summary_path) else NULL
test_text <- if (is.null(test_summary)) {
  "The test summary was not supplied."
} else {
  sprintf("%d tests; %d expectations passed; %d failures, %d errors, %d warnings, %d skips.",
          sum(test_summary$tests), sum(test_summary$passed), sum(test_summary$failed),
          sum(test_summary$errors), sum(test_summary$warnings), sum(test_summary$skipped))
}
package_check_path <- file.path(results_directory, "R-CMD-check.log")
package_check_lines <- if (file.exists(package_check_path)) readLines(package_check_path) else character()
package_check_status <- grep("^Status:", package_check_lines, value = TRUE)
package_check_text <- if (length(package_check_status)) tail(package_check_status, 1L) else
  "Status unavailable; package-check log was not supplied."
exact_directory <- if (dir.exists(arguments$exact)) normalizePath(arguments$exact) else NULL
exact_metrics <- if (!is.null(exact_directory) &&
                     file.exists(file.path(exact_directory, "metrics.csv"))) {
  utils::read.csv(file.path(exact_directory, "metrics.csv"), stringsAsFactors = FALSE)
} else NULL
exact_rejections <- if (!is.null(exact_metrics)) {
  utils::read.csv(file.path(exact_directory, "exact-rejection-rates.csv"))
} else NULL
edge_path <- file.path(results_directory, "edge-probes.csv")
edge_probes <- if (file.exists(edge_path)) utils::read.csv(edge_path, stringsAsFactors = FALSE) else NULL
if (is.null(edge_probes) ||
    !all(c("case", "family", "expected", "observed", "passed") %in% names(edge_probes)) ||
    !is.logical(edge_probes$passed) || anyNA(edge_probes$passed)) {
  stop("A complete edge-probes.csv is required before rendering the final report.",
       call. = FALSE)
}
edge_failures <- sum(!edge_probes$passed)
expected_error <- startsWith(edge_probes$expected, "error:")
separation_probe <- expected_error & edge_probes$family == "binomial" &
  grepl("^(complete_|quasi_)", edge_probes$case)
gaussian_exact_probe <- expected_error & edge_probes$family == "gaussian" &
  startsWith(edge_probes$case, "Gaussian_exact")
guard_group <- function(label, selected) {
  data.frame(Check = label, Cases = sum(selected),
             `Expected behavior met` = sum(edge_probes$passed[selected]),
             `Remaining failures` = sum(!edge_probes$passed[selected]),
             check.names = FALSE)
}
guard_table <- rbind(
  guard_group("Bernoulli separation examples", separation_probe),
  guard_group("Nonconstant exact Gaussian fits", gaussian_exact_probe),
  guard_group("Endpoint, constant-response, and other guards",
              expected_error & !separation_probe & !gaussian_exact_probe),
  guard_group("Interior/noisy positive controls", !expected_error)
)
comparison_path <- file.path(results_directory, "guard-comparison.csv")
guard_comparison <- if (file.exists(comparison_path)) {
  utils::read.csv(comparison_path, stringsAsFactors = FALSE)
} else NULL
comparison_section <- if (!is.null(guard_comparison)) {
  complete_match <- guard_comparison$value[
    guard_comparison$check == "identical_complete_replicate_table"
  ]
  compared <- guard_comparison$value[guard_comparison$check == "replicates_compared"]
  paste0('<p>An optional archived before/after comparison records ',
    if (length(complete_match) == 1L && complete_match == 1) 'identical complete replicate tables' else
      'changes that require inspection',
    if (length(compared) == 1L) paste0(' across ', format(compared, big.mark = ","), ' simulated datasets') else '',
    ' after adding the guards. ',
    '<a href="guard-comparison.csv">Read the recorded comparison</a>. ',
    'This comparison uses an archived earlier run; it is separate from the current-run reproduction commands below.</p>')
} else ''
failures <- replicates[replicates$status == "error", ]
failure_counts <- if (nrow(failures)) {
  stats::aggregate(list(datasets = failures$replicate),
                   by = list(case_id = failures$case_id, error = failures$error), FUN = length)
} else data.frame(case_id = character(), error = character(), datasets = integer())
null_table <- data.frame(
  Scenario = unname(scenario_labels[nulls$scenario]), n = nulls$n,
  `Rejected / requested` = paste0(nulls$rejections, " / ", nulls$requested),
  Errors = nulls$errors, `Rejection rate` = percentage(nulls$unconditional_rejection_rate),
  `95% pointwise CI` = interval_text(nulls$rejection_ci_lower, nulls$rejection_ci_upper),
  `Adjusted excess p` = format_p(nulls$bonferroni_excess_p_value),
  `Excess flag` = ifelse(nulls$excess_flag, "FLAG", "No"), check.names = FALSE
)
power_table <- data.frame(
  Family = power$family, n = power$n, gamma = power$gamma,
  `Rejected / requested` = paste0(power$rejections, " / ", power$requested),
  Errors = power$errors, Power = percentage(power$unconditional_rejection_rate),
  `95% pointwise CI` = interval_text(power$rejection_ci_lower, power$rejection_ci_upper),
  check.names = FALSE
)
stress_table <- data.frame(
  Scenario = unname(scenario_labels[stress$scenario]), n = stress$n,
  Requested = stress$requested, Successful = stress$successful, Errors = stress$errors,
  `Error rate` = percentage(stress$errors / stress$requested),
  `Unconditional rejection` = percentage(stress$unconditional_rejection_rate),
  `95% pointwise CI` = interval_text(stress$rejection_ci_lower, stress$rejection_ci_upper),
  `Rejection among successful calls` = percentage(stress$success_conditional_rejection_rate),
  check.names = FALSE
)
warnings_table <- data.frame(
  `Case ID` = summary$case_id, `Datasets with warnings` = summary$warning_datasets,
  check.names = FALSE
)
warnings_table <- warnings_table[warnings_table[[2L]] > 0L, , drop = FALSE]
format_run_time <- function(time_string) {
  format(as.POSIXct(time_string, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
         format = "%B %d, %Y at %I:%M:%S %p %Z", tz = "America/New_York")
}
flagged <- sum(nulls$excess_flag)
null_conclusion <- if (!flagged) {
  "None of the 12 regular null cells shows detected excess rejection after the Bonferroni adjustment."
} else {
  sprintf("%d of 12 regular null cells show detected excess rejection after the Bonferroni adjustment; inspect the flagged rows below.", flagged)
}
expected_constant_probability <- 0.98^40 + 0.02^40
regular_failures <- sum(summary$errors[summary$category != "stress"])
header_status <- if (edge_failures > 0L) {
  "Recorded edge-probe failures require review"
} else if (regular_failures == 0L && flagged == 0L) {
  "Recorded guard cases pass; Gaussian/Poisson study completed"
} else "Validation completed with findings requiring review"
exact_section <- if (is.null(exact_metrics)) {
  "<p>No exact-reference output was supplied.</p>"
} else {
  paste0("<p>The independent balanced-label check exhausts all 70 allocations of four treated labels among eight fixed outcomes. ",
    "Group mean differences provide the reference ranking. Public observed scores were checked for every allocation; ",
    "the full-orbit p-value helper was checked for all three alternatives. Nine seeded public Monte Carlo runs ",
    "(199 permutations each) were also compared draw by draw, including p-values and RNG state.</p>",
    html_table(exact_metrics), "<p>Exact rejection fractions over this fixed orbit:</p>",
    html_table(exact_rejections),
    '<p>References: <a href="../exact-results/metrics.csv">metrics</a>, ',
    '<a href="../exact-results/allocations.csv">all allocations</a>, ',
    '<a href="../exact-results/full-orbit-p-values.csv">full-orbit p-values</a>, ',
    '<a href="../exact-results/monte-carlo-runs.csv">Monte Carlo runs</a>.</p>')
}
edge_section <- paste0(
  "<p>Deterministic public-API tests cover boundary responses (all-zero counts, all-zero/all-one Bernoulli outcomes), ",
  "constant outcomes with constant offsets, row alignment, nuisance aliases, all supported families, fixed offsets, ",
  "and decreasing links. Constant outcomes with a varying offset can remain informative and are checked against ",
  "an independent fitted-model reference.</p>",
  if (!is.null(edge_probes)) paste0(
    '<p><a href="edge-probes.csv">Download all deterministic edge probes</a> · ',
    '<a href="edge-probes-session-info.txt">Probe session information</a></p>',
    html_table(edge_probes[, intersect(c("case", "family", "link", "expected",
      "observed", "passed", "p_value", "notes"), names(edge_probes)), drop = FALSE])) else
    "<p>Additional edge-probe details were not supplied to this report.</p>"
)
source_table <- data.frame(File = names(config$source_md5), MD5 = unname(config$source_md5),
                           check.names = FALSE)
reproduction <- sprintf(
  "R_LIBS=%s Rscript validation/validate-perglm.R --repetitions=%d --B=%d --workers=%d --output=validation/results\nR_LIBS=%s Rscript validation/exact-reference.R\nR_LIBS=%s Rscript validation/edge-probes.R\nRscript validation/summarize-validation.R --results=validation/results --exact=validation/exact-results",
  dirname(config$installed_package_path), config$requested_options$repetitions,
  config$requested_options$B, config$effective_workers,
  dirname(config$installed_package_path), dirname(config$installed_package_path)
)
guard_section <- paste0(
  '<section', if (edge_failures > 0L) ' class="unresolved"' else '',
  '><h2>Guard regression checks</h2><p><strong>',
  sum(edge_probes$passed), ' of ', nrow(edge_probes),
  ' recorded deterministic probes meet their stated expected behavior; ',
  edge_failures, ' remain failing.</strong> Expected errors and finite positive controls are both retained.</p>',
  '<p>Earlier validation identified complete/quasi separation in Bernoulli models and roundoff-only scores from nonconstant exact Gaussian fits. ',
  'Those defects motivated explicit fit-validation guards. The table below reports their current rerun status rather than carrying forward the earlier failure count.</p>',
  '<p>Fit-validation checks run once during null-score preparation. The permutation engine reuses the prepared quantities across all draws.</p>',
  '<ul><li><strong>Bernoulli separation:</strong> <code>check_glm_separation()</code> calls ',
  '<code>detectseparation::detect_separation()</code> using a logit geometry proxy and the column-normalized free design. ',
  'Fixed offsets are excluded from separation directions; the requested family, link, and offset remain in the actual model fit. ',
  'The result records <code>null_diagnostics$separation_checked</code>.</li>',
  '<li><strong>Gaussian residuals:</strong> <code>gaussian_fit_diagnostics()</code> compares fitted-response residual variation ',
  'with a bound derived from prediction-arithmetic roundoff, including nonconstant outcomes and fixed offsets. ',
  'The relative arithmetic envelope uses <code>gamma_k = k*u/(1-k*u)</code>, with <code>k = 4*(p+2)</code>, ',
  '<code>u = .Machine$double.eps/2</code>, and nuisance-design dimension <code>p</code>. ',
  'The result records <code>null_diagnostics$gaussian_residual_check</code>, including ',
  '<code>residual_norm</code>, <code>roundoff_bound</code>, and <code>relative_residual_norm</code>.</li></ul>',
  html_table(guard_table),
  if (edge_failures > 0L) paste0(
    '<p>The following recorded cases still fail their expected behavior. Favorable regular-case simulations do not resolve these failures.</p>',
    html_table(edge_probes[!edge_probes$passed,
      intersect(c("case", "family", "link", "expected", "observed", "p_value", "notes"),
                names(edge_probes)), drop = FALSE])) else
    '<p>All currently recorded guard and positive-control probes pass. This verifies the listed cases and does not establish coverage of every possible design or numerical boundary.</p>',
  comparison_section,
  '</section>'
)
html <- paste0('<!doctype html><html lang="en"><head><meta charset="utf-8">',
  '<meta name="viewport" content="width=device-width,initial-scale=1">',
  '<title>perglm_test validation</title><style>',
  'body{font-family:system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;color:#172b3a;background:#f2f6fa;margin:0;line-height:1.55}',
  'main{max-width:1180px;margin:auto;padding:36px 28px 64px}h1{font-size:34px;line-height:1.2;margin-bottom:10px}h2{margin-top:36px;font-size:23px}',
  'p{max-width:1050px}.subtitle{color:#4f6475}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:12px;margin:24px 0}',
  '.card,section{background:white;border:1px solid #dce5ed;border-radius:10px;padding:20px}.card strong{display:block;font-size:28px;color:#155981}.card span{color:#526879}',
  'section.unresolved{border-color:#d5a4a4;border-left:5px solid #a33131;background:#fff9f9}.unresolved h2{color:#922b2b}',
  'section{margin-top:18px}section h2{margin-top:0}.notice{border-left:4px solid #26749a;padding:12px 18px;background:#e9f3f8}',
  '.table-wrap{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:13px;margin:15px 0}th,td{text-align:left;border-bottom:1px solid #e4ebf1;padding:10px 9px;vertical-align:top}',
  'th{background:#eef4f8;color:#24465e;font-weight:650}tr:nth-child(even){background:#fafcfe}img{width:100%;height:auto;background:white}a{color:#155f8c}',
  'pre{overflow-x:auto;background:#102e43;color:#f2f7fb;padding:18px;border-radius:7px;font-size:12px}code{font-family:ui-monospace,SFMono-Regular,Consolas,monospace}',
  'details{margin:18px 0}summary{cursor:pointer;font-weight:600}footer{font-size:12px;color:#5f7484;margin-top:28px}',
  '@media(max-width:650px){main{padding:24px 14px}h1{font-size:28px}section{padding:15px}}</style></head><body><main>',
  '<h1><code>perglm_test()</code> validation</h1><p class="subtitle">', escape_html(header_status),
  '</p><p class="subtitle">Run completed ', escape_html(format_run_time(config$completed_utc)),
  '. Package ', escape_html(config$perglm_version), '; ', escape_html(config$R_version), '.</p>',
  '<div class="cards"><div class="card"><strong>', nrow(summary), '</strong><span>simulation cells</span></div>',
  '<div class="card"><strong>', format(sum(summary$requested), big.mark = ","), '</strong><span>requested datasets</span></div>',
  '<div class="card"><strong>', config$requested_options$B, '</strong><span>permutations per successful call</span></div>',
  '<div class="card"><strong>', edge_failures, '</strong><span>unresolved edge-probe failures</span></div></div>',
  '<p class="notice">This is a finite simulation study and an implementation audit. It does not prove calibration for all distributions, sample sizes, or model choices. ',
  'Every requested dataset is retained; failures give no decision and count as no rejection in unconditional rates.</p>',
  guard_section,
  '<section><h2>Implementation and exact-reference checks</h2><p>', escape_html(test_text), '</p>',
  '<p><code>R CMD check --no-manual --no-build-vignettes</code>: ', escape_html(package_check_text),
  ' <a href="R-CMD-check.log">Read the package-check log</a>.</p>',
  '<p>Existing public-API tests exercise default and custom fitting controls, batch-size and statistic-retention reproducibility, ',
  'and agreement of final RNG states with independent calculations. The simulation below uses the study-specific settings listed explicitly.</p>',
  '<p>Saved raw results, counts, rejection indicators, plus-one p-values, replicate seeds, and exact binomial intervals were rechecked. ',
  'All source-file hashes match the simulation snapshot.</p>', exact_section, '</section>',
  '<section><h2>Simulation design</h2>',
  '<p>Each cell uses ', config$requested_options$repetitions, ' independent datasets, a two-sided test, nominal level 0.05, and ',
  config$requested_options$B, ' independently sampled permutations. ',
  'Study-specific fitting controls are <code>epsilon = 1e-10</code> and <code>maxit = 100</code>, ',
  'with <code>alternative = "two.sided"</code>, <code>batch_size = 64</code>, and <code>keep_statistics = FALSE</code>. ',
  'This uses a tighter convergence tolerance and a larger iteration limit than the public defaults; the study evaluates this specified configuration. ',
  'There were ', regular_failures, ' errors across the ',
  format(sum(summary$requested[summary$category != "stress"]), big.mark = ","),
  ' regular-null/power datasets; the ', sum(stress$errors),
  ' sparse-exposure failures are reported separately below. ',
  'Two bounded nuisance covariates U and V are independent Uniform[-1,1]. ',
  'Ordinary exposures have prevalence 0.3 under marginal independence (MI), or Bernoulli probability ',
  'expit(logit(0.3) + 1.2U - 0.8V) in confounded conditional-mean (CM) cells. The fitted null is always <code>y ~ u + v</code>.</p>',
  '<p>Gaussian means are 0.4 + 0.6U - 0.4V + gamma X; Poisson log means are 0.2 + 0.5U - 0.3V + gamma X. ',
  'MI cases add an omitted 0.6U^2 term while keeping X independent of (Y,U,V). CM cases retain the correct fitted null mean. ',
  'The Gaussian heteroskedastic cell uses SD = exp(0.5U + 0.3X): its conditional mean null holds, while conditional independence is false because variance depends on X. ',
  'The Poisson variance-misspecification cell generates negative-binomial outcomes with size 2 and the correct Poisson working mean.</p>',
  '<p>Gaussian power coefficients are 0.25 and 0.5; Poisson coefficients are 0.3 and 0.6. Sample sizes are 100 and 500. ',
  'Separate sparse-exposure stress cases use n = 40 and independent X ~ Bernoulli(0.02).</p></section>',
  '<section><h2>Null calibration</h2><p>', escape_html(null_conclusion),
  ' The familywise check uses a one-sided binomial excess tail multiplied by 12; only adjusted p-values below 0.05 flag excess rejection. ',
  'The plotted intervals are 95% pointwise exact binomial intervals and are not simultaneous bands.</p>',
  embedded_image(null_plot, "Null rejection probabilities with exact pointwise confidence intervals"),
  '<p><a href="null-calibration.png">Download null-calibration figure</a> · ',
  '<a href="null-calibration-checks.csv">Download multiplicity checks</a></p>', html_table(null_table), '</section>',
  '<section><h2>Power</h2><p>The zero-effect points are the correctly specified confounded-CM null cells. ',
  'Each positive-effect cell has its own independent datasets and seeds. Intervals describe repeated-dataset rejection probabilities.</p>',
  embedded_image(power_plot, "Power versus exposure coefficient for two sample sizes"),
  '<p><a href="power.png">Download power figure</a></p>', html_table(power_table), '</section>',
  '<section><h2>Sparse-exposure failures and difficult inputs</h2>',
  '<p>With n = 40 and exposure prevalence 0.02, the probability that X is constant is ',
  percentage(expected_constant_probability), ' (0.98^40 + 0.02^40). Such datasets have no identifiable exposure contrast. ',
  'Observed errors and messages are reported below. Successful-call rejection rates are descriptive and do not establish validity conditional on successful completion.</p>',
  html_table(stress_table), html_table(failure_counts), edge_section,
  if (nrow(warnings_table)) paste0('<h3>Warnings retained</h3>', html_table(warnings_table)) else
    '<p>No simulation dataset produced a warning.</p>', '</section>',
  '<section><h2>Reproducibility and data</h2><p>RNG: ', escape_html(paste(config$rng_kind, collapse = " / ")),
  '. Base seed: ', config$base_seed, '. Data and permutation seeds are stored for every replicate. ',
  'The simulation uses the recorded installed package; reinstall this source version into that library before reproduction.</p>',
  '<pre>', escape_html(reproduction), '</pre>',
  '<p><a href="summary.csv">All 22 cell summaries</a> · <a href="replicates.csv.gz">Every replicate (compressed CSV)</a> · ',
  '<a href="scenarios.csv">Scenario specifications</a> · <a href="test-summary.csv">Automated-test summary</a> · ',
  '<a href="session-info.txt">R session information</a> · <a href="config.rds">Run configuration</a></p>',
  '<details><summary>Verified source hashes</summary>', html_table(source_table), '</details>',
  '<p>Installed simulation package: <code>', escape_html(config$installed_package_path), '</code>. ',
  'Simulation wall time: ', sprintf("%.1f", config$elapsed_seconds), ' seconds. Analyzer MD5: <code>',
  escape_html(unname(tools::md5sum(script_path))), '</code>.</p></section>',
  '<footer>Generated ', escape_html(format(Sys.time(), "%B %d, %Y at %I:%M:%S %p %Z", tz = "America/New_York")),
  '. Figures use base R graphics; no external assets are needed to view this report.</footer>',
  '</main></body></html>')
report_path <- file.path(results_directory, "report.html")
writeLines(html, report_path, useBytes = TRUE)
cat("Validated source snapshot and every raw/summary consistency check.\n")
cat(sprintf("Regular null excess flags after Bonferroni: %d / %d\n", flagged, nrow(nulls)))
cat(sprintf("Errors: %d regular/power; %d sparse-exposure stress.\n",
            regular_failures, sum(stress$errors)))
cat(sprintf("Unresolved deterministic edge-probe failures: %d.\n", edge_failures))
cat("Report: ", report_path, "\n", sep = "")
