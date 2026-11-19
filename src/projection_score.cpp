#include "studentize.h"

#include <R_ext/BLAS.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

const char* projection_error =
  "projection must be a valid object returned by prepare_projection().";
const char* cache_error =
  "projection contains invalid or incompatible cached factors.";
const char* exposure_error =
  "exposure must be a finite real numeric vector with one value per design row.";
const char* exposures_error =
  "exposures must be a finite real numeric matrix with one row per design row.";
const char* residual_error =
  "score_residuals must be a finite real numeric vector with one value per design row.";
const char* contribution_error =
  "Residualized score contributions must be finite.";
const char* alias_error =
  "exposure has no numerically identifiable variation after nuisance adjustment.";

void check_exposure_tolerance(const double tolerance, const Rcpp::NumericMatrix& z) {
  if (!std::isfinite(tolerance) || tolerance < 0 || tolerance >= 1) {
    Rcpp::stop("exposure_tol must be finite, at least zero and less than one.");
  }
  if (tolerance > 0) {
    // Centering an exposure is valid only when an intercept is in the design.
    for (int j = 0; j < z.ncol(); ++j) {
      bool intercept = true;
      for (int i = 0; i < z.nrow(); ++i) {
        if (z(i, j) != 1) {
          intercept = false;
          break;
        }
      }
      if (intercept) {
        return;
      }
    }
    Rcpp::stop("An intercept is required when exposure_tol is positive.");
  }
}

bool numeric_matrix(SEXP value) {
  return perglm::is_real_numeric(value) && Rf_isMatrix(value);
}

bool numeric_vector(SEXP value, R_xlen_t length) {
  return perglm::is_real_numeric(value) &&
    Rf_isNull(Rf_getAttrib(value, R_DimSymbol)) && XLENGTH(value) == length;
}

void require_finite(const double* values, R_xlen_t length,
                    const char* message) {
  for (R_xlen_t i = 0; i < length; ++i) {
    if (i > 0 && (i & 65535) == 0) {
      Rcpp::checkUserInterrupt();
    }
    if (!std::isfinite(values[i])) {
      Rcpp::stop(message);
    }
  }
}

// The factors are checked once per call, then reused for every exposure column.
// Coercing integer matrices/vectors also translates integer NA to real NA.
class Projection {
 public:
  explicit Projection(SEXP value) : z(0, 0), q(0, 0), r(0, 0),
                                    sqrt_weights(0), pivot(0), n(0), p(0) {
    if (!Rf_inherits(value, "perglm_projection") || TYPEOF(value) != VECSXP) {
      Rcpp::stop(projection_error);
    }
    Rcpp::List cache(value);
    if (!cache.containsElementNamed("z") ||
        !cache.containsElementNamed("q") ||
        !cache.containsElementNamed("r") ||
        !cache.containsElementNamed("sqrt_weights") ||
        !cache.containsElementNamed("qr")) {
      Rcpp::stop(cache_error);
    }
    SEXP z_value = cache["z"];
    SEXP q_value = cache["q"];
    SEXP r_value = cache["r"];
    SEXP weight_value = cache["sqrt_weights"];
    SEXP qr_value = cache["qr"];
    if (!numeric_matrix(z_value) || !numeric_matrix(q_value) ||
        !numeric_matrix(r_value) || TYPEOF(qr_value) != VECSXP) {
      Rcpp::stop(cache_error);
    }
    z = Rcpp::as<Rcpp::NumericMatrix>(z_value);
    q = Rcpp::as<Rcpp::NumericMatrix>(q_value);
    r = Rcpp::as<Rcpp::NumericMatrix>(r_value);
    n = z.nrow();
    p = z.ncol();
    if (n < 2 || p < 1 || p > n || q.nrow() != n || q.ncol() != p ||
        r.nrow() != p || r.ncol() != p || !numeric_vector(weight_value, n)) {
      Rcpp::stop(cache_error);
    }
    sqrt_weights = Rcpp::as<Rcpp::NumericVector>(weight_value);
    require_finite(z.begin(), z.size(), cache_error);
    require_finite(q.begin(), q.size(), cache_error);
    require_finite(r.begin(), r.size(), cache_error);
    require_finite(sqrt_weights.begin(), n, cache_error);
    for (int i = 0; i < n; ++i) {
      if (sqrt_weights[i] <= 0) {
        Rcpp::stop(cache_error);
      }
    }
    for (int j = 0; j < p; ++j) {
      // QR diagonal entries may legitimately be negative.
      if (r(j, j) == 0) {
        Rcpp::stop(cache_error);
      }
    }

    Rcpp::List qr(qr_value);
    if (!qr.containsElementNamed("pivot") ||
        !qr.containsElementNamed("rank")) {
      Rcpp::stop(cache_error);
    }
    SEXP pivot_value = qr["pivot"];
    SEXP rank_value = qr["rank"];
    if (TYPEOF(pivot_value) != INTSXP || XLENGTH(pivot_value) != p ||
        TYPEOF(rank_value) != INTSXP || XLENGTH(rank_value) != 1 ||
        INTEGER(rank_value)[0] != p) {
      Rcpp::stop(cache_error);
    }
    pivot = Rcpp::IntegerVector(pivot_value);
    std::vector<bool> seen(p, false);
    for (int j = 0; j < p; ++j) {
      const int index = pivot[j];
      if (index == NA_INTEGER || index < 1 || index > p || seen[index - 1]) {
        Rcpp::stop(cache_error);
      }
      seen[index - 1] = true;
    }
  }

  Rcpp::NumericMatrix z;
  Rcpp::NumericMatrix q;
  Rcpp::NumericMatrix r;
  Rcpp::NumericVector sqrt_weights;
  Rcpp::IntegerVector pivot;
  int n;
  int p;
};

class ScoreWorkspace {
 public:
  explicit ScoreWorkspace(const Projection& projection) :
    weighted(projection.n), coefficients(projection.p),
    beta(projection.p), contributions(projection.n) {}

  double compute(const Projection& projection, const double* exposure,
                 const double* score_residuals, const char* invalid_error,
                 const double exposure_tol) {
    const int n = projection.n;
    const int p = projection.p;
    double minimum = exposure[0];
    double maximum = exposure[0];
    for (int i = 0; i < n; ++i) {
      if (i > 0 && (i & 65535) == 0) {
        Rcpp::checkUserInterrupt();
      }
      if (!std::isfinite(exposure[i])) {
        Rcpp::stop(invalid_error);
      }
      minimum = std::min(minimum, exposure[i]);
      maximum = std::max(maximum, exposure[i]);
    }
    const double range = maximum - minimum;
    const double scale = std::max(std::abs(minimum), std::abs(maximum));
    if (exposure_tol > 0 && minimum == maximum) {
      Rcpp::stop(alias_error);
    }
    for (int i = 0; i < n; ++i) {
      double value = exposure[i];
      if (exposure_tol > 0) {
        // A positive affine transformation leaves the studentized score
        // unchanged with an intercept. Subtract before scaling to retain
        // variation in exposures with a large common offset. Min/max are
        // unchanged by permutations, so the same labels get the same values.
        value = std::isfinite(range) ? (value - minimum) / range - 0.5 :
          (value / scale - minimum / scale) /
            (maximum / scale - minimum / scale) - 0.5;
      }
      weighted[i] = projection.sqrt_weights[i] * value;
      contributions[i] = value;
      if (!std::isfinite(weighted[i])) {
        Rcpp::stop("Exposure residualization produced nonfinite values.");
      }
    }

    const double one = 1;
    const double zero = 0;
    const double negative_one = -1;
    const int increment = 1;
    const double exposure_norm = exposure_tol > 0 ?
      F77_CALL(dnrm2)(&n, weighted.data(), &increment) : 0;
    // Q and R describe sqrt(W) Z[, pivot]. No factorization is repeated.
    F77_CALL(dgemv)("T", &n, &p, &one, projection.q.begin(), &n,
                    weighted.data(), &increment, &zero, coefficients.data(),
                    &increment FCONE);
    F77_CALL(dtrsv)("U", "N", "N", &p, projection.r.begin(), &p,
                    coefficients.data(), &increment FCONE FCONE FCONE);
    for (int j = 0; j < p; ++j) {
      if (!std::isfinite(coefficients[j])) {
        Rcpp::stop("Exposure residualization produced nonfinite values.");
      }
      beta[projection.pivot[j] - 1] = coefficients[j];
    }
    F77_CALL(dgemv)("N", &n, &p, &negative_one, projection.z.begin(), &n,
                    beta.data(), &increment, &one, contributions.data(),
                    &increment FCONE);
    for (int i = 0; i < n; ++i) {
      if (!std::isfinite(contributions[i])) {
        Rcpp::stop("Exposure residualization produced nonfinite values.");
      }
      if (exposure_tol > 0) {
        weighted[i] = projection.sqrt_weights[i] * contributions[i];
      }
    }
    if (exposure_tol > 0 &&
        (exposure_norm == 0 || F77_CALL(dnrm2)(&n, weighted.data(), &increment) <=
          exposure_tol * exposure_norm)) {
      Rcpp::stop(alias_error);
    }
    for (int i = 0; i < n; ++i) {
      contributions[i] *= score_residuals[i];
    }
    return perglm::studentize_column(contributions.data(), n,
                                     contribution_error);
  }

 private:
  std::vector<double> weighted;
  std::vector<double> coefficients;
  std::vector<double> beta;
  std::vector<double> contributions;
};

}  // namespace

// Residualize one exposure and studentize its product with fixed score residuals.
// Positive exposure_tol enables GLM identifiability checks and normalization;
// the default retains the general projection primitive without centering.
// [[Rcpp::export(rng = false)]]
double residualized_score_cpp(SEXP projection, SEXP exposure,
                             SEXP score_residuals, double exposure_tol = 0) {
  Projection factors(projection);
  check_exposure_tolerance(exposure_tol, factors.z);
  if (!numeric_vector(exposure, factors.n)) {
    Rcpp::stop(exposure_error);
  }
  if (!numeric_vector(score_residuals, factors.n)) {
    Rcpp::stop(residual_error);
  }
  Rcpp::Shield<SEXP> numeric_exposure(Rf_coerceVector(exposure, REALSXP));
  Rcpp::Shield<SEXP> numeric_residuals(Rf_coerceVector(score_residuals, REALSXP));
  require_finite(REAL(numeric_residuals), factors.n, residual_error);
  ScoreWorkspace workspace(factors);
  return workspace.compute(factors, REAL(numeric_exposure),
                            REAL(numeric_residuals), exposure_error, exposure_tol);
}

// Process a block using fixed factors/residuals and reusable column workspaces.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector residualized_score_batch_cpp(SEXP projection,
                                                SEXP exposures,
                                                SEXP score_residuals,
                                                double exposure_tol = 0) {
  Projection factors(projection);
  check_exposure_tolerance(exposure_tol, factors.z);
  if (!numeric_matrix(exposures)) {
    Rcpp::stop(exposures_error);
  }
  const SEXP dimensions = Rf_getAttrib(exposures, R_DimSymbol);
  const R_xlen_t n = INTEGER(dimensions)[0];
  const R_xlen_t columns = INTEGER(dimensions)[1];
  if (n != factors.n) {
    Rcpp::stop(exposures_error);
  }
  if (!numeric_vector(score_residuals, n)) {
    Rcpp::stop(residual_error);
  }
  Rcpp::Shield<SEXP> numeric_exposures(Rf_coerceVector(exposures, REALSXP));
  Rcpp::Shield<SEXP> numeric_residuals(Rf_coerceVector(score_residuals, REALSXP));
  require_finite(REAL(numeric_residuals), n, residual_error);
  ScoreWorkspace workspace(factors);
  Rcpp::NumericVector result(columns);
  const double* values = REAL(numeric_exposures);
  for (R_xlen_t column = 0; column < columns; ++column) {
    if ((column & 63) == 0) {
      Rcpp::checkUserInterrupt();
    }
    try {
      result[column] = workspace.compute(factors, values + column * n,
                                         REAL(numeric_residuals), exposures_error,
                                         exposure_tol);
    } catch (const Rcpp::exception& error) {
      Rcpp::stop("Column %lld: %s", static_cast<long long>(column + 1), error.what());
    }
  }
  const SEXP dimnames = Rf_getAttrib(exposures, R_DimNamesSymbol);
  if (!Rf_isNull(dimnames)) {
    const SEXP column_names = VECTOR_ELT(dimnames, 1);
    if (!Rf_isNull(column_names)) {
      result.attr("names") = column_names;
    }
  }
  return result;
}
