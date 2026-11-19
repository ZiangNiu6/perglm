# Run from the package root against its installed development version:
# R_LIBS=/path/to/library Rscript validation/exact-reference.R
# This diagnostic does not use perglm's projection or studentization helpers.
# It conditions on fixed outcomes and exhausts the balanced binary label orbit.

library(perglm)

p_value_helper <- getFromNamespace("permutation_p_value", "perglm")
output_directory <- file.path("validation", "exact-results")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

n <- 8L
treated <- 4L
assignments <- utils::combn(n, treated)
orbit_size <- ncol(assignments)
exposures <- matrix(0, nrow = n, ncol = orbit_size)
for (column in seq_len(orbit_size)) {
  exposures[assignments[, column], column] <- 1
}
keys <- apply(exposures, 2L, paste0, collapse = "")

# Dyadic values avoid incidental floating ties: every four-element subset has
# a unique sum. Centering preserves that property and gives fixed real outcomes.
y <- 2^(seq_len(n) - 1L) / 64
y <- y - mean(y)
data <- data.frame(y = y, row.names = paste0("observation_", seq_len(n)))
centered_y <- y - mean(y)

# With four labels of each kind, residual exposure is X - 1/2. Thus
# S = (n/4) * (mean(Y_1) - mean(Y_0)), Q = sum((Y-mean(Y))^2)/4,
# and T = S / sqrt(Q - S^2/n). This includes the empirical n divisor.
differences <- vapply(seq_len(orbit_size), function(column) {
  x <- exposures[, column]
  mean(y[x == 1]) - mean(y[x == 0])
}, numeric(1))
numerators <- (n / 4) * differences
uncentered_squares <- sum(centered_y^2) / 4
reference_statistics <- numerators /
  sqrt(uncentered_squares - numerators^2 / n)

public_statistics <- numeric(orbit_size)
public_B1_p_values <- numeric(orbit_size)
observed_seed <- 7241L
for (column in seq_len(orbit_size)) {
  set.seed(observed_seed)
  result <- perglm_test(y ~ 1, data = data, exposure = exposures[, column],
                        family = stats::gaussian(), B = 1L)
  public_statistics[column] <- unname(result$statistic)
  public_B1_p_values[column] <- result$p.value
}

allocation_results <- data.frame(
  allocation = seq_len(orbit_size), labels = keys,
  treated_positions = apply(assignments, 2L, paste, collapse = ";"),
  group_mean_difference = differences,
  reference_statistic = reference_statistics,
  public_statistic = public_statistics,
  absolute_error = abs(public_statistics - reference_statistics),
  public_B1_p_value = public_B1_p_values
)

alternatives <- c("greater", "less", "two.sided")
tail_values <- function(values, alternative) {
  if (alternative == "two.sided") abs(values) else {
    if (alternative == "less") -values else values
  }
}

orbit_results <- do.call(rbind, lapply(alternatives, function(alternative) {
  group_scores <- tail_values(differences, alternative)
  group_p <- vapply(seq_len(orbit_size), function(column) {
    sum(group_scores >= group_scores[column]) / orbit_size
  }, numeric(1))
  helper_reference_p <- vapply(seq_len(orbit_size), function(column) {
    p_value_helper(reference_statistics[column], reference_statistics[-column],
                    alternative = alternative)
  }, numeric(1))
  helper_public_p <- vapply(seq_len(orbit_size), function(column) {
    p_value_helper(public_statistics[column], public_statistics[-column],
                    alternative = alternative)
  }, numeric(1))
  data.frame(allocation = seq_len(orbit_size), labels = keys,
               alternative = alternative, group_p_value = group_p,
               reference_score_p_value = helper_reference_p,
               public_score_p_value = helper_public_p,
               reference_p_matches = helper_reference_p == group_p,
               public_p_matches = helper_public_p == group_p)
}))

alpha_grid <- c(0.05, 0.1, 0.2)
rejection_results <- do.call(rbind, lapply(alternatives, function(alternative) {
  results <- orbit_results[orbit_results$alternative == alternative, ]
  do.call(rbind, lapply(alpha_grid, function(alpha) {
    group_rate <- mean(results$group_p_value <= alpha)
    reference_rate <- mean(results$reference_score_p_value <= alpha)
    public_rate <- mean(results$public_score_p_value <= alpha)
    data.frame(alternative = alternative, alpha = alpha,
                 group_rejection_fraction = group_rate,
                 reference_rejection_fraction = reference_rate,
                 public_rejection_fraction = public_rate,
                 all_at_most_alpha = max(group_rate, reference_rate, public_rate) <=
                   alpha + 1e-14)
  }))
}))

mc_budget <- 199L
observed_allocations <- c(1L, 23L, 70L)
mc_seeds <- c(119L, 223L, 337L)
mc_results <- list()
draw_results <- list()
run <- 0L
for (index in seq_along(observed_allocations)) {
  allocation <- observed_allocations[index]
  x <- exposures[, allocation]
  seed <- mc_seeds[index]
  set.seed(seed)
  sampled_keys <- vapply(seq_len(mc_budget), function(draw) {
    paste0(x[sample.int(n, n, replace = FALSE)], collapse = "")
  }, character(1))
  sampled_allocations <- match(sampled_keys, keys)
  reference_rng <- .Random.seed
  expected_statistics <- reference_statistics[sampled_allocations]
  expected_differences <- differences[sampled_allocations]

  for (alternative in alternatives) {
    run <- run + 1L
    set.seed(seed)
    result <- perglm_test(y ~ 1, data = data, exposure = x,
                          family = stats::gaussian(), B = mc_budget,
                          alternative = alternative, batch_size = 17L,
                          keep_statistics = TRUE)
    actual_rng <- .Random.seed
    group_observed <- tail_values(differences[allocation], alternative)
    group_permuted <- tail_values(expected_differences, alternative)
    expected_exceedances <- sum(group_permuted >= group_observed)
    expected_p <- (1 + expected_exceedances) / (mc_budget + 1)
    reference_p <- p_value_helper(reference_statistics[allocation],
                                   expected_statistics, alternative)
    actual_statistics <- unname(result$permuted_statistics)
    mc_results[[run]] <- data.frame(
      observed_allocation = allocation, labels = keys[allocation], seed = seed,
      alternative = alternative, B = mc_budget,
      reference_observed_statistic = reference_statistics[allocation],
      public_observed_statistic = unname(result$statistic),
      expected_exceedances = expected_exceedances,
      public_exceedances = result$exceedances,
      group_p_value = expected_p, reference_score_p_value = reference_p,
      public_p_value = result$p.value,
      p_values_match = expected_p == reference_p && expected_p == result$p.value,
      rng_states_match = identical(actual_rng, reference_rng),
      all_assignments_mapped = all(!is.na(sampled_allocations)),
      max_absolute_statistic_error = max(abs(actual_statistics - expected_statistics)),
      max_scaled_statistic_error = max(abs(actual_statistics - expected_statistics) /
                                        pmax(1, abs(expected_statistics)))
    )
    draw_results[[run]] <- data.frame(
      observed_allocation = allocation, seed = seed, alternative = alternative,
      draw = seq_len(mc_budget), sampled_allocation = sampled_allocations,
      labels = sampled_keys, reference_statistic = expected_statistics,
      public_statistic = actual_statistics,
      absolute_error = abs(actual_statistics - expected_statistics)
    )
  }
}
mc_results <- do.call(rbind, mc_results)
draw_results <- do.call(rbind, draw_results)

tolerance <- 1e-11
maximum_scaled_observed_error <- max(
  abs(public_statistics - reference_statistics) /
    pmax(1, abs(reference_statistics))
)
metrics <- data.frame(
  check = c("balanced_allocations", "distinct_signed_group_differences",
              "distinct_absolute_group_differences", "observed_score_scaled_error",
              "reference_orbit_p_mismatches", "public_orbit_p_mismatches",
              "exact_rejection_fraction_violations", "mc_score_scaled_error",
              "mc_p_value_mismatches", "mc_rng_mismatches", "mc_unmapped_runs"),
  value = c(orbit_size, length(unique(differences)), length(unique(abs(differences))),
              maximum_scaled_observed_error, sum(!orbit_results$reference_p_matches),
              sum(!orbit_results$public_p_matches),
              sum(!rejection_results$all_at_most_alpha),
              max(mc_results$max_scaled_statistic_error),
              sum(!mc_results$p_values_match), sum(!mc_results$rng_states_match),
              sum(!mc_results$all_assignments_mapped)),
  passed = c(orbit_size == 70L, length(unique(differences)) == 70L,
               length(unique(abs(differences))) == 35L,
               maximum_scaled_observed_error <= tolerance,
               all(orbit_results$reference_p_matches), all(orbit_results$public_p_matches),
               all(rejection_results$all_at_most_alpha),
               all(mc_results$max_scaled_statistic_error <= tolerance),
               all(mc_results$p_values_match), all(mc_results$rng_states_match),
               all(mc_results$all_assignments_mapped))
)

utils::write.csv(data.frame(observation = rownames(data), y = y),
                 file.path(output_directory, "fixed-outcomes.csv"), row.names = FALSE)
utils::write.csv(allocation_results, file.path(output_directory, "allocations.csv"),
                 row.names = FALSE)
utils::write.csv(orbit_results, file.path(output_directory, "full-orbit-p-values.csv"),
                 row.names = FALSE)
utils::write.csv(rejection_results, file.path(output_directory, "exact-rejection-rates.csv"),
                 row.names = FALSE)
utils::write.csv(mc_results, file.path(output_directory, "monte-carlo-runs.csv"),
                 row.names = FALSE)
utils::write.csv(draw_results, file.path(output_directory, "monte-carlo-draws.csv"),
                 row.names = FALSE)
utils::write.csv(metrics, file.path(output_directory, "metrics.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_directory, "session-info.txt"))

print(metrics, row.names = FALSE)
print(rejection_results, row.names = FALSE)
print(mc_results, row.names = FALSE)
if (!all(metrics$passed)) {
  stop("Exact-reference diagnostics failed; inspect the saved CSV comparisons.",
         call. = FALSE)
}
