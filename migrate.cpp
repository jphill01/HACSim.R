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
    double m;
    IntegerVector ind1, ind2;
    
    if (m != 0) {
        if (model == "Step") {
            for (i = 1; i < (K - 1); i++) {
                for (j = 1; j < (K - 1); j++) {
                        ind1 = RcppArmadillo::sample(perms, ceil(perms * m / 2), false);
                        ind2 = RcppArmadillo::sample(perms, ceil(perms * m / 2), false);
                        tmp = pop[ind1, RcppArmadillo::sample(i, ind1, true)];
                        pop[ind1, RcppArmadillo::sample(i, ind1, true)] = pop[ind2, RcppArmadillo::sample(j, ind2, true)];
                        pop[ind2, RcppArmadillo::sample(j, ind2, true)] = tmp;

                    }
                }
            }
        }
    return pop;
}
