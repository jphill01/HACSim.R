// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_OPENMP_WARNING
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <set>
using namespace Rcpp;

// [[Rcpp::export]]
arma::Cube<int> migrate(arma::Cube<int> pop) {
    
    String model;
    double m = 0;
    int N = 0, K = 0;
    
    if (m != 0 && model == "Step") {
        for (int i = 0; i < (K - 1); i++) {
            IntegerVector ind = RcppArmadillo::sample<IntegerVector>(pop[i], ceil(N * m), true);
            int tmp = pop[ind[i]];
            pop[ind[i]] = pop[ind[i + 1]];
            pop[ind[i + 1]] = tmp;
        }
        return pop;
        } else {
            return pop;
            }
    }
