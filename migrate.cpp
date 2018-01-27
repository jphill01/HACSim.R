// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_OPENMP_WARNING
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <set>
using namespace Rcpp;

double m;
String model;
int K, num_specs, perms, tmp;

IntegerVector ind;

// [[Rcpp::export]]
arma::Cube<int> migrate(arma::Cube<int> pop) {
        if (model == "Island") {
            ind = RcppArmadillo::sample<IntegerVector>(perms, ceil(num_specs * m), true);
                for (int i = 0; i < (K - 1); i++) {
                    tmp = pop[ind, i];
                    pop[ind, i] = pop[ind, i + 1];
                    pop[ind, i + 1] = tmp;
                }
            }
            return pop;
} else {
    return pop
}
