// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_OPENMP_WARNING
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <set>
using namespace Rcpp;

// [[Rcpp::export]]
arma::Cube<int> migrate(arma::Cube<int> pop) {
    
    // String model;
    double m = 0;
    int K = 0, perms = 0;
    // N = 0;
    // int num_specs = ceil(N / K);
    
    if (m == 0) { // K == 1 && model == NULL
        return pop;
    } else {
        IntegerVector ind = RcppArmadillo::sample<IntegerVector>(perms, ceil(perms * m), false);
        for (int i = 0; i < (K - 1); i++) {     // K > 1 && m == 0 && model == "Step"
            int tmp = pop[ind[i]];
            pop[ind[i]] = pop[ind[i + 1]];
            pop[ind[i + 1]] = tmp;
        }
        return pop;
    }
    
    }

