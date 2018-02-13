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
    int K = 0, perms;
    int N;
    int num_specs = ceil(N / K);
    
    if (K = 1, m = 0 && model == NULL) {
        return pop;
    } else {
        for (int i = 0; i < (K - 1); i++) {     // K > 1 && m != 0 && model == "Step"
            IntegerVector ind = RcppArmadillo::sample<IntegerVector>(perms, ceil(num_specs * m), false);
            int tmp = pop[ind[i]];
            pop[ind[i]] = pop[ind[i + 1]];
            pop[ind[i + 1]] = tmp;
        }
        return pop;
    }
    
    }
