// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_OPENMP_WARNING
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <set>
using namespace Rcpp;

double m;
String model;
int K, perms, tmp;

IntegerVector ind;
arma::Cube<int> pop;

// [[Rcpp::export]]
arma::Cube<int> migrate() {
    if (m != 0) {
        if (model == "Step") {
            ind = RcppArmadillo::sample<IntegerVector>(perms, ceil(perms * m/2), false);
            for (int i = 0; i < (K - 2); i++) {
                        tmp = pop[ind, i];
                        pop[ind, i] = pop[ind, i + 1];
                        pop[ind, i + 1] = tmp;
                    }
                }
            }
    return pop;
}
