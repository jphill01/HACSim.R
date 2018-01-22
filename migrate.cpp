// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_OPENMP_WARNING
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <set>
using namespace Rcpp;

// [[Rcpp::export]]
const arma::Cube<int> pop;

arma::Cube<int> migrate(pop) {
    int i, j, k, K, m;
    std::string model;
    IntegerVector ind1, ind2;
    
    if (m != 0) {
        if (model == "Step") {
            for (i = 0; i < K; i++) {
                for (j = 0; j < K; j++) {
                    for (k = 1; k < (K - 1); k++) {
                        ind1 = sample(perms, ceiling(perms * m / 2), false);
                        ind2 = sample(perms, ceiling(perms * m / 2), false);
                        tmp = pop[ind1, RcppArmadillo::sample(k)];
                        pop[ind1, RcppArmadillo::sample(k)] = pop[ind2, RcppArmadillo::sample(k)];
                        pop[ind2, RcppArmadillo::sample(k)] = tmp;
                        return pop;
                    }
                }
            }
        }
    }
}
