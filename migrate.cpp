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
    int K = 0, N = 0;
    int num_specs = ceil(N / K);
    IntegerVector haps;
    NumericVector probs;
    
    if (m == 0) { // K == 1 && model == NULL
        return pop;
    } else {
        IntegerVector ind = RcppArmadillo::sample<IntegerVector>(haps, ceil(num_specs * m), false);
        for (int i = 0; i < (K - 1); i++) {     // K > 1 && m == 0 && model == "Step"
            int tmp = pop[ind[i]];
            pop[ind[i]] = pop[ind[i + 1]];
            pop[ind[i + 1]] = tmp;
        }

    // Update probabilities to account for migrated individuals
        if (K > 1 && m != 0) {
            probs = (1 - m) * probs + m * probs[ind]
            probs = probs / sum(probs)
        }
        return pop;
    }
    
    }

