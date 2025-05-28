#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::export]]
DataFrame Metro_Hast(double lambda_start, double rho_start, double t_start, int iter, NumericVector a, NumericVector x) {
  NumericVector lambda(iter, NA_REAL);
  NumericVector rho(iter, NA_REAL);
  NumericVector t(iter, NA_REAL);
  NumericVector loglik(iter, NA_REAL);
  IntegerVector accept(iter, NA_INTEGER);
  
  lambda[0] = lambda_start;
  rho[0] = rho_start;
  t[0] = t_start;
  
  if (lambda[0] < 0 || lambda[0] > 0.1 || rho[0] < 1 || rho[0] > 15 || t[0] < 5 || t[0] > 65) {
    loglik[0] = -INFINITY;
  } else {
    NumericVector liknum(a.size());
    for (int j = 0; j < a.size(); j++) {
      if (t[0] < a[j]) {
        liknum[j] = log(exp((lambda[0] * (t[0] - a[j]) - t[0] * rho[0] * lambda[0]) * fabs(x[j] - 1)) * 
          pow(1 - exp(lambda[0] * (t[0] - a[j]) - t[0] * rho[0] * lambda[0]), x[j]));
      } else {
        liknum[j] = log(exp(-rho[0] * lambda[0] * a[j] * fabs(x[j] - 1)) * 
          pow(1 - exp(-rho[0] * lambda[0] * a[j]), x[j]));
      }
    }
    loglik[0] = sum(liknum);
  }
  accept[0] = 0;
  
  NumericVector prop_sd = {0.0001, 0.5, 10};
  
  for (int i = 1; i < iter; i++) {
    NumericVector mu = {lambda[i - 1], rho[i - 1], t[i - 1]};
    NumericVector draws = {1,1,1};
    NumericVector theta_new = { R::rnorm(mu[0],prop_sd[0]),  R::rnorm(mu[1],prop_sd[1]), R::rnorm(mu[2],prop_sd[2]) };
    
    
    double loglik_new;
    if (theta_new[0] < 0 || theta_new[0] > 0.1 || theta_new[1] < 1 || theta_new[1] > 15 || theta_new[2] < 5 || theta_new[2] > 65) {
      loglik_new = -INFINITY;
    } else {
      NumericVector liknum(a.size());
      for (int j = 0; j < a.size(); j++) {
        if (theta_new[2] < a[j]) {
          liknum[j] = log(exp((theta_new[0] * (theta_new[2] - a[j]) - theta_new[2] * theta_new[1] * theta_new[0]) * fabs(x[j] - 1)) * 
            pow(1 - exp(theta_new[0] * (theta_new[2] - a[j]) - theta_new[2] * theta_new[1] * theta_new[0]), x[j]));
        } else {
          liknum[j] = log(exp(-theta_new[1] * theta_new[0] * a[j] * fabs(x[j] - 1)) * 
            pow(1 - exp(-theta_new[1] * theta_new[0] * a[j]), x[j]));
        }
      }
      loglik_new = sum(liknum);
    }
    
    double log_ratio = loglik_new - loglik[i - 1];
    double log_accept_prob = std::min(log_ratio, 0.0);
    
    if (log(R::runif(0, 1)) < log_accept_prob) {
      lambda[i] = theta_new[0];
      rho[i] = theta_new[1];
      t[i] = theta_new[2];
      loglik[i] = loglik_new;
      accept[i] = 1;
    } else {
      lambda[i] = lambda[i - 1];
      rho[i] = rho[i - 1];
      t[i] = t[i - 1];
      loglik[i] = loglik[i - 1];
      accept[i] = 0;
    }
  }
  return DataFrame::create(_["lambda"] = lambda, _["rho"] = rho, _["t"] = t, _["loglik"] = loglik, _["accept"] = accept);
}
