#include <Rcpp.h>
#include "cpp_helpers.h"  
using namespace Rcpp;

// Forward declare internal helpers
double calc_objective_sum_c(const Rcpp::NumericMatrix &dist, const Rcpp::IntegerVector &selected, double lambda_var);
double calc_objective_maxmin_c(const Rcpp::NumericMatrix &dist, const Rcpp::IntegerVector &selected);

// ---------- Exported wrappers ----------

// [[Rcpp::export]]
double calc_objective_sum(Rcpp::NumericMatrix dist, Rcpp::IntegerVector selected, double lambda_var = 0.0) {
    return calc_objective_sum_c(dist, selected, lambda_var);
}

// [[Rcpp::export]]
double calc_objective_maxmin(Rcpp::NumericMatrix dist, Rcpp::IntegerVector selected) {
    return calc_objective_maxmin_c(dist, selected);
}

// [[Rcpp::export]]
Rcpp::List local_search_swap(Rcpp::NumericMatrix dist, Rcpp::IntegerVector selected, 
                             Rcpp::IntegerVector candidates, std::string objective, 
                             int max_iter, double lambda_var = 0.0) {
    int n_sel = selected.size();
    int n_cand = candidates.size();
    double current_obj = -R_PosInf;

    if(objective == "sum") {
        current_obj = calc_objective_sum_c(dist, selected, lambda_var);
    } else {
        current_obj = calc_objective_maxmin_c(dist, selected);
    }

    bool improved = true;
    int iter = 0;

    while(improved && iter < max_iter) {
        improved = false;
        for(int i = 0; i < n_sel; ++i) {
            for(int j = 0; j < n_cand; ++j) {
                double new_obj;
                Rcpp::IntegerVector temp_selected = clone(selected);
                temp_selected[i] = candidates[j];

                if(objective == "sum") {
                    new_obj = calc_objective_sum_c(dist, temp_selected, lambda_var);
                } else {
                    new_obj = calc_objective_maxmin_c(dist, temp_selected);
                }

                if(new_obj > current_obj) {
                    int temp = selected[i];
                    selected[i] = candidates[j];
                    candidates[j] = temp;
                    current_obj = new_obj;
                    improved = true;
                    break;
                }
            }
            if(improved) break;
        }
        ++iter;
    }

    return Rcpp::List::create(
        Rcpp::Named("selected") = selected,
        Rcpp::Named("objective") = current_obj,
        Rcpp::Named("iterations") = iter
    );
}

