.glm_null_data <- function(family_name = "gaussian", seed = 431L, n = 150L) {
  set.seed(seed)
  z1 <- seq(-1, 1, length.out = n)
  z2 <- cos(seq_len(n) * 0.37)
  group <- factor(rep(c("first", "second", "third"), length.out = n),
                  levels = c("first", "second", "third"))
  linear <- 0.2 + 0.25 * z1 - 0.15 * z2 +
    0.1 * (group == "second") - 0.1 * (group == "third")
  mu <- exp(linear)
  y <- switch(family_name,
    gaussian = linear + rnorm(n, sd = 0.6),
    binomial = rbinom(n, 1, plogis(linear)),
    poisson = rpois(n, mu),
    Gamma = rgamma(n, shape = 5, scale = mu / 5),
    inverse.gaussian = {
      shape <- 5
      v <- rnorm(n)^2
      candidate <- mu + mu^2 * v / (2 * shape) -
        mu / (2 * shape) * sqrt(4 * mu * shape * v + mu^2 * v^2)
      ifelse(runif(n) <= mu / (mu + candidate), candidate, mu^2 / candidate)
    },
    negative.binomial = rnbinom(n, mu = mu, size = 1.7)
  )
  data.frame(y, z1, z2, group,
              row.names = sprintf("sample_%03d", seq_len(n)))
}

test_that("the wrapper matches independently fitted standard null models", {
  families <- list(gaussian = gaussian(), binomial = binomial(),
                    poisson = poisson(), Gamma = Gamma("log"),
                    inverse.gaussian = inverse.gaussian("log"))
  for (name in names(families)) {
    data <- .glm_null_data(name)
    family <- families[[name]]
    control <- glm.control(epsilon = 1e-10, maxit = 100)
    reference_fit <- glm(y ~ z1 + z2 + group, data = data, family = family,
                         control = control, x = TRUE, y = TRUE, model = TRUE,
                         na.action = na.fail, singular.ok = FALSE)
    reference <- prepare_glm_score(reference_fit)
    prepared <- fit_glm_null(y ~ z1 + z2 + group, data = data,
                             family = family, control = control)
    expect_s3_class(prepared, "perglm_glm_score")
    expect_identical(prepared$family_name, name)
    expect_equal(coef(prepared$fit), coef(reference_fit), tolerance = 1e-12)
    expect_equal(prepared$fit$fitted.values, reference_fit$fitted.values,
                 tolerance = 1e-12)
    expect_equal(prepared$score_residuals, reference$score_residuals,
                 tolerance = 1e-12)
    expect_equal(prepared$working_weights, reference$working_weights,
                 tolerance = 1e-12)
    expect_identical(prepared$row_names, rownames(data))
    expect_identical(rownames(prepared$fit$x), rownames(data))
    expect_identical(names(prepared$fit$y), rownames(data))
    expect_equal(unname(prepared$fit$prior.weights), rep(1, nrow(data)))
    expect_true(is.data.frame(prepared$fit$model))
    expect_identical(prepared$fit$model$group, data$group)
    expect_true(all(c("groupsecond", "groupthird") %in% colnames(prepared$fit$x)))
  }
})

test_that("family objects, stock constructors and exact family names agree", {
  constructors <- list(gaussian = stats::gaussian, binomial = stats::binomial,
                        poisson = stats::poisson, Gamma = stats::Gamma,
                        inverse.gaussian = stats::inverse.gaussian)
  for (name in names(constructors)) {
    data <- .glm_null_data(name)
    constructor <- constructors[[name]]
    from_name <- fit_glm_null(y ~ z1 + z2 + group, data, family = name)
    from_function <- fit_glm_null(y ~ z1 + z2 + group, data, family = constructor)
    from_object <- fit_glm_null(y ~ z1 + z2 + group, data, family = constructor())
    expect_equal(coef(from_name$fit), coef(from_object$fit), tolerance = 1e-12)
    expect_equal(coef(from_function$fit), coef(from_object$fit), tolerance = 1e-12)
    expect_equal(from_name$score_residuals, from_object$score_residuals,
                 tolerance = 1e-12)
    expect_equal(from_function$working_weights, from_object$working_weights,
                 tolerance = 1e-12)
  }
  data <- .glm_null_data()
  default <- fit_glm_null(y ~ z1 + z2, data)
  explicit <- fit_glm_null(y ~ z1 + z2, data, family = gaussian())
  expect_identical(default$family_name, "gaussian")
  expect_equal(coef(default$fit), coef(explicit$fit), tolerance = 1e-12)

  poisson_data <- .glm_null_data("poisson")
  factory <- function(...) stats::poisson(...)
  from_factory <- fit_glm_null(y ~ z1 + z2, poisson_data, family = factory)
  from_stock <- fit_glm_null(y ~ z1 + z2, poisson_data, family = poisson())
  expect_equal(coef(from_factory$fit), coef(from_stock$fit), tolerance = 1e-12)
  expect_equal(from_factory$score_residuals, from_stock$score_residuals,
               tolerance = 1e-12)
})

test_that("estimated and fixed negative-binomial sizes follow distinct fitting paths", {
  skip_if_not_installed("MASS")
  data <- .glm_null_data("negative.binomial", n = 180L)
  control <- glm.control(epsilon = 1e-10, maxit = 100)
  reference_fit <- MASS::glm.nb(y ~ z1 + z2 + group, data = data,
                               control = control, x = TRUE, y = TRUE, model = TRUE,
                               na.action = na.fail)
  reference <- prepare_glm_score(reference_fit)
  estimated <- fit_glm_null(y ~ z1 + z2 + group, data,
                            family = "negative.binomial", control = control)
  expect_identical(class(estimated$fit), c("negbin", "glm", "lm"))
  expect_identical(estimated$family_name, "negative.binomial")
  expect_equal(estimated$theta, reference_fit$theta, tolerance = 1e-12)
  expect_equal(coef(estimated$fit), coef(reference_fit), tolerance = 1e-12)
  expect_equal(estimated$score_residuals, reference$score_residuals,
               tolerance = 1e-12)
  expect_false("family" %in% names(as.list(estimated$fit$call)))
  updated <- update(estimated$fit, . ~ . - group)
  update_reference <- MASS::glm.nb(y ~ z1 + z2, data = data,
                                    control = control, x = TRUE, y = TRUE,
                                    model = TRUE, na.action = na.fail)
  # update.negbin() warm-starts size, so independent fits can differ slightly.
  expect_equal(coef(updated), coef(update_reference), tolerance = 1e-9)
  expect_equal(updated$theta, update_reference$theta, tolerance = 1e-9)

  fixed_family <- MASS::negative.binomial(theta = 1.7, link = "log")
  fixed_reference <- glm(y ~ z1 + z2 + group, data = data, family = fixed_family,
                         control = control, x = TRUE, y = TRUE, model = TRUE,
                         na.action = na.fail, singular.ok = FALSE)
  fixed <- fit_glm_null(y ~ z1 + z2 + group, data,
                        family = fixed_family, control = control)
  expect_identical(class(fixed$fit), c("glm", "lm"))
  expect_equal(fixed$theta, 1.7)
  expect_equal(coef(fixed$fit), coef(fixed_reference), tolerance = 1e-12)
  expect_equal(fixed$working_weights,
               prepare_glm_score(fixed_reference)$working_weights,
               tolerance = 1e-12)
})

test_that("partial controls are normalized and stored calls remain reusable", {
  analysis_data <- .glm_null_data("poisson")
  null_formula <- y ~ z1 + z2
  working_family <- poisson()
  fitting_control <- list(epsilon = 1e-10, maxit = 100)
  prepared <- fit_glm_null(null_formula, data = analysis_data,
                           family = working_family, control = fitting_control)
  expect_identical(prepared$fit$control,
                    do.call(stats::glm.control, fitting_control))
  expect_identical(prepared$fit$call$formula, quote(null_formula))
  expect_identical(prepared$fit$call$data, quote(analysis_data))
  expect_identical(prepared$fit$call$family, quote(working_family))
  expect_true("fitting_control" %in% all.vars(prepared$fit$call$control))
  expect_identical(eval(prepared$fit$call$control),
                    do.call(stats::glm.control, fitting_control))
  updated <- update(prepared$fit, . ~ . + group)
  reference <- glm(y ~ z1 + z2 + group, data = analysis_data,
                   family = working_family,
                   control = do.call(stats::glm.control, fitting_control),
                   x = TRUE, y = TRUE, model = TRUE)
  expect_equal(coef(updated), coef(reference), tolerance = 1e-12)
  expect_s3_class(prepare_glm_score(updated), "perglm_glm_score")
})

test_that("fixed formula offsets and data-frame subclasses are supported", {
  data <- .glm_null_data("poisson")
  data$fixed_offset <- 0.1 * sin(seq_len(nrow(data)))
  class(data) <- c("analysis_frame", class(data))
  original <- data
  prepared <- fit_glm_null(y ~ z1 + z2 + group + offset(fixed_offset),
                           data, family = poisson())
  reference <- glm(y ~ z1 + z2 + group + offset(fixed_offset), data = data,
                   family = poisson(), x = TRUE, y = TRUE, model = TRUE)
  expect_equal(coef(prepared$fit), coef(reference), tolerance = 1e-12)
  expect_equal(prepared$fit$offset, reference$offset, tolerance = 1e-12)
  expect_identical(prepared$row_names, rownames(data))
  expect_identical(data, original)
})

test_that("missing data cannot silently change the inference rows", {
  data <- .glm_null_data("poisson")
  for (field in c("y", "z1", "group")) {
    incomplete <- data
    incomplete[[field]][1L] <- NA
    expect_error(fit_glm_null(y ~ z1 + z2 + group, incomplete, family = poisson()))
  }
  expect_error(fit_glm_null(y ~ missing_predictor, data, family = poisson()))
  expect_error(fit_glm_null(missing_response ~ z1, data, family = poisson()))
})

test_that("the wrapper rejects invalid formulas, data and model ranks", {
  data <- .glm_null_data()
  for (formula in list(NULL, 1, "y ~ z1", ~ z1, y ~ z1 ~ z2)) {
    expect_error(fit_glm_null(formula, data))
  }
  for (invalid_data in list(NULL, as.matrix(data), as.list(data),
                            data[FALSE, ], data[1L, , drop = FALSE])) {
    expect_error(fit_glm_null(y ~ z1 + z2, invalid_data))
  }
  expect_error(fit_glm_null(y ~ z1))
  expect_error(fit_glm_null(y ~ z1 + z2 - 1, data))
  expect_error(fit_glm_null(y ~ z1 + z2 + I(2 * z1), data))
})

test_that("unsupported families and fitting controls are rejected", {
  data <- .glm_null_data("poisson")
  for (family in list(NULL, "normal", "Poisson", "nb", "quasipoisson",
                      "negative.binomial()", c("gaussian", "poisson"),
                      list(), quasipoisson())) {
    expect_error(fit_glm_null(y ~ z1 + z2, data, family = family))
  }
  custom_link <- poisson()
  custom_link$link <- "custom_log"
  expect_error(fit_glm_null(y ~ z1 + z2, data, family = custom_link))
  for (control in list(NULL, 1, "default", list(unknown_control = TRUE),
                       list(epsilon = 0), list(maxit = 0), list(epsilon = NA_real_))) {
    expect_error(fit_glm_null(y ~ z1 + z2, data, family = poisson(), control = control))
  }
  for (tol in list(0, 1, NA_real_, Inf, c(0.1, 0.2), 1i)) {
    expect_error(fit_glm_null(y ~ z1 + z2, data, family = poisson(), tol = tol))
  }
})

test_that("the wrapper does not accept row-changing or weighted extra arguments", {
  data <- .glm_null_data("poisson")
  expect_error(fit_glm_null(y ~ z1 + z2, data, weights = rep(1, nrow(data))),
               "unused argument")
  expect_error(fit_glm_null(y ~ z1 + z2, data, subset = seq_len(100)),
               "unused argument")
  expect_error(fit_glm_null(y ~ z1 + z2, data, na.action = na.omit),
               "unused argument")
  expect_error(fit_glm_null(y ~ z1 + z2, data, x = FALSE), "unused argument")
  expect_error(fit_glm_null(y ~ z1 + z2, data, y = FALSE), "unused argument")
  expect_error(fit_glm_null(y ~ z1 + z2, data, model = FALSE), "unused argument")
})

test_that("a fit stopped before convergence cannot produce a prepared null", {
  data <- .glm_null_data("poisson")
  expect_error(suppressWarnings(fit_glm_null(y ~ z1 + z2 + group, data,
                                             family = poisson(),
                                             control = glm.control(maxit = 1))))
})

test_that("fitted null quantities feed the native score kernels", {
  if (!exists("_perglm_residualized_score_cpp", inherits = TRUE)) {
    skip("Native registration is unavailable in source-only adapter checks.")
  }
  data <- .glm_null_data("poisson")
  prepared <- fit_glm_null(y ~ z1 + z2 + group, data, family = poisson())
  exposures <- cbind(binary = as.numeric(sin(seq_len(nrow(data)) * sqrt(2)) > 0),
                     continuous = sin(seq_len(nrow(data))))
  adjusted_norms <- apply(exposures, 2L, function(exposure) {
    sqrt(sum(residualize_exposure(prepared$projection, exposure)^2))
  })
  expect_true(all(adjusted_norms > 0.1))
  reference <- apply(exposures, 2L, function(exposure) {
    studentize_score(residualize_exposure(prepared$projection, exposure) *
                       prepared$score_residuals)
  })
  batch <- residualized_score_batch_cpp(prepared$projection, exposures,
                                        prepared$score_residuals)
  expect_equal(batch, reference, tolerance = 1e-11)
  expect_equal(residualized_score_cpp(prepared$projection, exposures[, 1L],
                                      prepared$score_residuals),
               unname(reference[1L]), tolerance = 1e-11)
})
