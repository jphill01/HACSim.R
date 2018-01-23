// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_OPENMP_WARNING
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <set>
using namespace Rcpp;

// [[Rcpp::export]]
arma::Cube<int> migrate(arma::Cube<int> pop) {
    String model;
    int i, j, k, K, perms, tmp;
    double m = 0;
    IntegerVector ind1, ind2;
    
    if (m != 0) {
        if (model == "Step") {
            for (i = 0; i < K; i++) {
                for (j = 0; j < K; j++) {
                    for (k = 1; k < (K - 1); k++) {
                        ind1 = sample(perms, ceil(perms * m / 2), false);
                        ind2 = sample(perms, ceil(perms * m / 2), false);
                        tmp = pop[ind1, RcppArmadillo::sample(k, k, true)];
                        pop[ind1, RcppArmadillo::sample(k, k, true)] = pop[ind2, RcppArmadillo::sample(k, k, true)];
                        pop[ind2, RcppArmadillo::sample(k, k, true)] = tmp;
                        }
                    }
                }
            }
        }
    return pop;
}
