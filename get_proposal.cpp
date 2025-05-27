#include <Rcpp.h>
using namespace Rcpp;

// Function to calculate the value of the log likelihood
// [[Rcpp::export]]
double calcLogLik(double lambda, double rho, double t, NumericVector a, IntegerVector x) {
  if (lambda < 0 || lambda > 0.1 || rho < 1 || rho > 15 || t < 5 || t > 65) {
    return R_NegInf;
  }
  
  int n = a.size();
  NumericVector lik(n);
  
  for (int i = 0; i < n; i++) {
    if (t < a[i]) {
      lik[i] = log(exp((lambda * (t - a[i]) - t * rho * lambda) * std::abs(x[i] - 1)) *
        pow(1 - exp(lambda * (t - a[i]) - t * rho * lambda), x[i]));
    } else {
      lik[i] = log(exp(-rho * lambda * a[i] * std::abs(x[i] - 1)) *
        pow(1 - exp(-rho * lambda * a[i]), x[i]));
    }
  }
  
  return sum(lik);
}

// [[Rcpp::export]]
double generate_rnorm(int n, double mean, double sd) {
  return Rcpp::rnorm(n, mean, sd);
}

// Function to run the Metropolis-Hastings Algorithm
// [[Rcpp::export]]
NumericMatrix get_proposal(double lambda_start, double rho_start, double t_start, int iter, NumericVector a, IntegerVector x) {
  NumericVector lambda(iter), rho(iter), t(iter), lik(iter), accept(iter);
  lambda[0] = lambda_start;
  rho[0] = rho_start;
  t[0] = t_start;
  lik[0] = calcLogLik(lambda_start, rho_start, t_start, a, x);
  accept[0] = 0;
  
  NumericVector prop_sd = {0.0001, 0.5, 10};
  NumericVector numdraws = {1,1,1};

  
  for (int i = 1; i < iter; i++) {
    NumericVector mu = {lambda[i - 1], rho[i - 1], t[i - 1]};
    NumericVector theta_new = generate_rnorm(numdraws,mu, prop_sd);
    
    double loglik_new = calcLogLik(theta_new[0], theta_new[1], theta_new[2], a, x);
    double loglik_old = calcLogLik(lambda[i - 1], rho[i - 1], t[i - 1], a, x);
    double log_ratio = loglik_new - loglik_old;
    double log_accept_prob = std::min(log_ratio, 0.0);
    
    if (log(R::runif(0, 1)) < log_accept_prob) {
      lambda[i] = theta_new[0];
      rho[i] = theta_new[1];
      t[i] = theta_new[2];
      lik[i] = loglik_new;
      accept[i] = 1;
    } else {
      lambda[i] = lambda[i - 1];
      rho[i] = rho[i - 1];
      t[i] = t[i - 1];
      lik[i] = loglik_old;
      accept[i] = 0;
    }
  }
  
  double p_cur = sum(accept) / iter;
  NumericMatrix prop_sd_new = cov(cbind(lambda, rho, t)) * R::qnorm(0.238 / 2, 0, 1) / R::qnorm(p_cur / 2);
  
  return prop_sd_new;
}
