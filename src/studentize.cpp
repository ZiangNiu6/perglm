#include "studentize.h"

#include <algorithm>
#include <cmath>

namespace perglm {

const char* vector_error =
  "contributions must be a finite real numeric vector of length at least two.";
const char* matrix_error =
  "contributions must be a finite real numeric matrix with at least two rows.";
const char* variance_error =
  "Score contributions must have positive empirical variance.";

// Preserve R's numeric contract for classed objects such as factors and dates.
// Ordinary vectors and matrices avoid a callback to R.
bool is_real_numeric(SEXP value) {
  if (TYPEOF(value) != REALSXP && TYPEOF(value) != INTSXP) {
    return false;
  }
  if (OBJECT(value)) {
    Rcpp::Function is_numeric = Rcpp::Environment::base_env()["is.numeric"];
    return Rcpp::as<bool>(is_numeric(value));
  }
  return true;
}

// Neumaier compensation helps even on platforms where long double and double
// have the same precision. No temporary score vectors are needed.
class CompensatedSum {
 public:
  CompensatedSum() : total_(0), correction_(0) {}

  void add(double value) {
    const long double term = value;
    const long double next = total_ + term;
    if (std::fabs(total_) >= std::fabs(term)) {
      correction_ += (total_ - next) + term;
    } else {
      correction_ += (term - next) + total_;
    }
    total_ = next;
  }

  long double value() const {
    return total_ + correction_;
  }

 private:
  long double total_;
  long double correction_;
};

double studentize_column(const double* values, R_xlen_t n,
                         const char* invalid_error) {
  double scale = 0;
  double difference_scale = 0;
  bool finite_differences = true;
  const double anchor = values[0];

  for (R_xlen_t i = 0; i < n; ++i) {
    if (i > 0 && (i & 65535) == 0) {
      Rcpp::checkUserInterrupt();
    }
    const double value = values[i];
    if (!std::isfinite(value)) {
      Rcpp::stop(invalid_error);
    }
    scale = std::max(scale, std::fabs(value));
    const double difference = value - anchor;
    if (std::isfinite(difference)) {
      difference_scale = std::max(difference_scale, std::fabs(difference));
    } else {
      finite_differences = false;
    }
  }
  if (scale == 0 || (finite_differences && difference_scale == 0)) {
    Rcpp::stop(variance_error);
  }

  CompensatedSum numerator;
  CompensatedSum centered_sum;
  for (R_xlen_t i = 0; i < n; ++i) {
    const double scaled = values[i] / scale;
    numerator.add(scaled);
    centered_sum.add(finite_differences ?
      (values[i] - anchor) / difference_scale : scaled);
  }
  const double center = static_cast<double>(centered_sum.value() / n);

  CompensatedSum sum_squares;
  for (R_xlen_t i = 0; i < n; ++i) {
    const double value = finite_differences ?
      (values[i] - anchor) / difference_scale : values[i] / scale;
    const double centered = value - center;
    sum_squares.add(centered * centered);
  }
  const long double ss = sum_squares.value();
  if (ss == 0) {
    Rcpp::stop(variance_error);
  }

  const double relative_scale = finite_differences ? difference_scale / scale : 1;
  return static_cast<double>(numerator.value() /
    (relative_scale * std::sqrt(ss)));
}

}  // namespace perglm

// Internal scalar counterpart of the R reference studentize_score().
// [[Rcpp::export(rng = false)]]
double studentize_score_cpp(SEXP contributions) {
  if (!perglm::is_real_numeric(contributions) ||
      !Rf_isNull(Rf_getAttrib(contributions, R_DimSymbol)) ||
      XLENGTH(contributions) < 2) {
    Rcpp::stop(perglm::vector_error);
  }

  // Coercion is a no-op for doubles and maps integer NA to real NA correctly.
  Rcpp::Shield<SEXP> numeric(Rf_coerceVector(contributions, REALSXP));
  return perglm::studentize_column(REAL(numeric), XLENGTH(numeric),
                                   perglm::vector_error);
}

// Each matrix column is an independent contribution vector. One R/C++ crossing
// processes the whole block, with no copies or allocations per column.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector studentize_score_batch_cpp(SEXP contributions) {
  if (!perglm::is_real_numeric(contributions) || !Rf_isMatrix(contributions)) {
    Rcpp::stop(perglm::matrix_error);
  }
  const SEXP dimensions = Rf_getAttrib(contributions, R_DimSymbol);
  const R_xlen_t n = INTEGER(dimensions)[0];
  const R_xlen_t columns = INTEGER(dimensions)[1];
  if (n < 2) {
    Rcpp::stop(perglm::matrix_error);
  }

  Rcpp::Shield<SEXP> numeric(Rf_coerceVector(contributions, REALSXP));
  const double* values = REAL(numeric);
  Rcpp::NumericVector result(columns);
  for (R_xlen_t column = 0; column < columns; ++column) {
    if ((column & 63) == 0) {
      Rcpp::checkUserInterrupt();
    }
    try {
      result[column] = perglm::studentize_column(values + column * n, n,
                                                 perglm::matrix_error);
    } catch (const Rcpp::exception& error) {
      Rcpp::stop("Column %lld: %s", static_cast<long long>(column + 1), error.what());
    }
  }

  const SEXP dimnames = Rf_getAttrib(contributions, R_DimNamesSymbol);
  if (!Rf_isNull(dimnames)) {
    const SEXP column_names = VECTOR_ELT(dimnames, 1);
    if (!Rf_isNull(column_names)) {
      result.attr("names") = column_names;
    }
  }
  return result;
}
