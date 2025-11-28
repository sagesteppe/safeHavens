// No Rcpp::export here; only internal helpers
#include <Rcpp.h>
#include "cpp_helpers.h" 
using namespace Rcpp;

// ---------- Helper non-exported functions ----------
double calc_objective_sum_c(const NumericMatrix &dist, const IntegerVector &selected, double lambda_var = 0.0) {
    int n = selected.size();
    if(n < 2) return 0.0;

    std::vector<double> dists;
    dists.reserve((n * (n - 1)) / 2);

    for(int i = 0; i < n - 1; ++i) {
        for(int j = i + 1; j < n; ++j) {
            int ii = selected[i] - 1;
            int jj = selected[j] - 1;
            dists.push_back(dist(ii, jj));
        }
    }

    double obj_disp = 0.0;
    for(size_t k = 0; k < dists.size(); ++k) obj_disp += dists[k];

    double obj_var = 0.0;
    if(dists.size() >= 2) {
        double mean = obj_disp / dists.size();
        double sum_sq_diff = 0.0;
        for(size_t k = 0; k < dists.size(); ++k) {
            double diff = dists[k] - mean;
            sum_sq_diff += diff * diff;
        }
        double variance = sum_sq_diff / (dists.size() - 1);
        obj_var = -variance;
    }

    return obj_disp + lambda_var * obj_var;
}

double calc_objective_maxmin_c(const NumericMatrix &dist, const IntegerVector &selected) {
    int n = selected.size();
    if(n < 2) return 0.0;

    double min_dist = R_PosInf;
    for(int i = 0; i < n - 1; ++i) {
        for(int j = i + 1; j < n; ++j) {
            int ii = selected[i] - 1;
            int jj = selected[j] - 1;
            double d = dist(ii, jj);
            if(d < min_dist) min_dist = d;
        }
    }
    return min_dist;
}

