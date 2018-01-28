// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_OPENMP_WARNING
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <set>
using namespace Rcpp;

// [[Rcpp::export]]
arma::Cube<int> migrate(arma::Cube<int> pop) {
    String model;
    double m;
    int K, perms, tmp;
    
    if (model == "Island") {
        IntegerVector ind = RcppArmadillo::sample<IntegerVector>(perms, ceil(perms * m), true);
        for (int i = 0; i < (K - 1); i++) {
            tmp = pop[ind, i];
            pop[ind, i] = pop[ind, i + 1];
            pop[ind, i + 1] = tmp;
        }
        return pop;
        } else {
            return pop;
            }
    }
