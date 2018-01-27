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
    int K, perms;
    
    if (model == "Island") {
        for (int i = 0; i < (K - 1); i++) {
            IntegerVector ind = RcppArmadillo::sample<IntegerVector>(perms, ceil(perms * m), true);
            int tmp = pop[ind, i];
            pop[ind, i] = pop[ind, i + 1];
            pop[ind, i + 1] = tmp;
        }
        return pop;
        } else {
            return pop;
            }
    }
