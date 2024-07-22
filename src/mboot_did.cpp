#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector mboot_did(NumericVector linrep, int nboot) {
  double k1 = 0.5 * (1 - sqrt(5.0));
  double k2 = 0.5 * (1 + sqrt(5.0));
  double pkappa = 0.5 * (1 + sqrt(5.0)) / sqrt(5.0);
  
  int n = linrep.size();
  NumericVector boot_did(nboot);

  for (int i = 0; i < nboot; ++i) {
    NumericVector v = rbinom(n, 1, pkappa);
    for (int j = 0; j < n; ++j) {
      v[j] = v[j] == 1 ? k1 : k2;
    } 
    boot_did[i] = mean(linrep * v);
  }

  return boot_did;
}

