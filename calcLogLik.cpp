#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
double calcLogLik(double lambda, double rho, double t, NumericVector a, NumericVector x) {
  if (lambda < 0 || lambda > 0.1 || rho < 1 || rho > 15 || t < 5 || t > 65) {
    return R_NegInf;
  }
  
  int n = a.size();
  NumericVector lik(n);
  
  for (int i = 0; i < n; i++) {
    if (t < a[i]) {
      lik[i] = std::log(std::exp((lambda * (t - a[i]) - t * rho * lambda) * std::abs(x[i] - 1)) *
        std::pow(1 - std::exp(lambda * (t - a[i]) - t * rho * lambda), x[i]));
    } else {
      lik[i] = std::log(std::exp(-rho * lambda * a[i] * std::abs(x[i] - 1)) *
        std::pow(1 - std::exp(-rho * lambda * a[i]), x[i]));
    }
  }
  
  return sum(lik);
}
