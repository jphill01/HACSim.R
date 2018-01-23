// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_OPENMP_WARNING
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <set>
using namespace Rcpp;

arma::Cube<int> pop;
double m;
const String model;
int i, j, k, K, perms, tmp;

// [[Rcpp::export]]
arma::Cube<int> migrate_cpp(arma::Cube<int> pop) {
    if (m != 0) {
        if (model == "Step") {
            for (i = 0; i < K; i++) {
                for (j = 0; j < K; j++) {
                    for(k = 1; k < (K - 1); k++) {
                        i = RcppArmadillo::sample(perms, ceil(perms * m / 2), false);
                        j = RcppArmadillo::sample(perms, ceil(perms * m / 2), false);
                        tmp = pop[i, RcppArmadillo::sample(k, K, true)];
                        pop[i, RcppArmadillo::sample(k, K, true)] = pop[j, RcppArmadillo::sample(k, K, true)];
                        pop[j, RcppArmadillo::sample(k, K, true)] = tmp;
                    }
                }
            }
        }
    }
    return pop;
}
