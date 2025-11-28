#ifndef CPP_HELPERS_H
#define CPP_HELPERS_H

#include <Rcpp.h>

// Declare helper functions
double calc_objective_sum_c(const Rcpp::NumericMatrix &dist, 
                           const Rcpp::IntegerVector &selected, 
                           double lambda_var);

double calc_objective_maxmin_c(const Rcpp::NumericMatrix &dist, 
                               const Rcpp::IntegerVector &selected);

#endif
