// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_OPENMP_WARNING
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <set>
using namespace Rcpp;

arma::Cube<int> pop;
double m;
const String model;
int i, ind, K, perms;
double tmp;

// [[Rcpp::export]]
arma::Cube<int> migrate_cpp(arma::Cube<int> pop) {
    if (m != 0) {
        if (model == "Step") {
            ind <- sample(perms, size = ceiling(perms * m/2), replace = FALSE)
            for (i = 0; i < (K - 1); i++) {
                        tmp = pop[ind, i];
                        pop[ind, i] = pop[ind,, i + 1];
                        pop[ind, i + 1] = tmp;
                    }
                }
            }
    return pop;
}
