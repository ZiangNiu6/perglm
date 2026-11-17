#ifndef PERGLM_STUDENTIZE_H
#define PERGLM_STUDENTIZE_H

#include <Rcpp.h>

namespace perglm {

bool is_real_numeric(SEXP value);
double studentize_column(const double* values, R_xlen_t n,
                         const char* invalid_error);

}  // namespace perglm

#endif
