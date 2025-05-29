
#include <cmath>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Function to compute log-likelihood
double compute_loglik(double lambda, double rho, double t, NumericVector a, NumericVector x) {
  if (lambda < 0 || lambda > 0.1 || rho < 1 || rho > 15 || t < 5 || t > 65) {
    return -INFINITY;
  }
  
  double loglik_sum = 0.0;
  for (int j = 0; j < a.size(); j++) {
    double liknum;
    if (t < a[j]) {
      liknum = log(exp((lambda * (t - a[j]) - t * rho * lambda) * std::abs(x[j] - 1)) *
                     pow(1 - exp(lambda * (t - a[j]) - t * rho * lambda), x[j]));
    } else {
      liknum = log(exp(-rho * lambda * a[j] * std::abs(x[j] - 1)) *
                     pow(1 - exp(-rho * lambda * a[j]), x[j]));
    }
    loglik_sum += liknum;
  }
  return loglik_sum;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector mvrnorm_rcpp(int n, arma::vec mu, arma::mat Sigma) {
  int d = mu.n_elem;
  arma::mat Y = arma::randn(d, n); // Generate standard normal variables
  arma::mat L = arma::chol(Sigma); // Cholesky decomposition
  arma::vec sample = mu + L * Y;   // Generate multivariate normal sample
  
  return Rcpp::wrap(sample); // Convert arma::vec to Rcpp::NumericVector
}

// [[Rcpp::export]]
NumericMatrix compute_cov_matrix(NumericVector lambda, NumericVector rho, NumericVector t, int burn_iter, double p_cur) {
  arma::mat data(burn_iter, 3);
  
  // Fill matrix with parameter values
  for (int i = 0; i < burn_iter; i++) {
    data(i, 0) = lambda[i];
    data(i, 1) = rho[i];
    data(i, 2) = t[i];
  }
  
  // Compute covariance matrix and scale by quantiles
  arma::mat proposal_cov = arma::cov(data) * (R::qnorm(0.238 / 2, 0, 1, 1, 0) / R::qnorm(p_cur / 2, 0, 1, 1, 0));
  
  // Convert arma::mat to Rcpp::NumericMatrix
  return Rcpp::wrap(proposal_cov);
}

// [[Rcpp::export]]

arma::mat convert_to_arma(NumericMatrix numeric_matrix) {
  arma::mat arma_mat(numeric_matrix.begin(), numeric_matrix.nrow(), numeric_matrix.ncol(), false);
  return arma_mat;
}

arma::vec convert_to_arma_vec(NumericVector numeric_vector) {
  arma::vec arma_vector(numeric_vector.begin(), numeric_vector.size(), false);
  return arma_vector;
}


// [[Rcpp::export]]
DataFrame Metro_Hast_cpp(double lambda_start, double rho_start, double t_start,
                     int burn_iter, int iter, NumericVector a, NumericVector x) {
  NumericVector lambda(burn_iter + iter, NA_REAL);
  NumericVector rho(burn_iter + iter, NA_REAL);
  NumericVector t(burn_iter + iter, NA_REAL);
  NumericVector loglik(burn_iter + iter, NA_REAL);
  NumericVector accept(burn_iter + iter, NA_REAL);
  
  lambda[0] = lambda_start;
  rho[0] = rho_start;
  t[0] = t_start;
  loglik[0] = compute_loglik(lambda[0], rho[0], t[0], a, x);
  accept[0] = 0;
  
  NumericVector prop_sd = {0.001, 0.5, 10};
  
  // Burn-in phase
  for (int i = 1; i < burn_iter; i++) {
    NumericVector theta_new = {R::rnorm(lambda[i - 1], prop_sd[0]),
      R::rnorm(rho[i - 1], prop_sd[1]),
      R::rnorm(t[i - 1], prop_sd[2])};
    
    double loglik_new = compute_loglik(theta_new[0], theta_new[1], theta_new[2], a, x);
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
  
  // Compute new proposal covariance matrix
  double p_cur = sum(accept[Rcpp::Range(0, burn_iter - 1)]) / burn_iter;
  NumericMatrix proposal_cov = compute_cov_matrix(lambda, rho, t, burn_iter, p_cur);
  arma::mat Sigma = convert_to_arma(proposal_cov);
  
  
  
  // Post burn-in phase using multivariate normal sampling
  for (int i = burn_iter; i < (burn_iter + iter); i++) {
    arma::vec mu = convert_to_arma_vec({lambda[i - 1], rho[i - 1], t[i - 1]});
    NumericVector theta_new = mvrnorm_rcpp(1, mu, Sigma);
    
    double loglik_new = compute_loglik(theta_new[0], theta_new[1], theta_new[2], a, x);
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
  
  return DataFrame::create(
    _["lambda"] = lambda,
    _["rho"] = rho,
    _["t"] = t,
    _["loglik"] = loglik,
    _["accept"] = accept
  );
}
