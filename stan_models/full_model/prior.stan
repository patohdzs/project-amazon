data {
  int<lower=0> S; // Number of sites
  int<lower=0> K_theta; // Number of coefficients on theta
  int<lower=0> K_gamma; // Number of coefficients on gamma
  int<lower=0> N_theta;
  int<lower=0> N_gamma;

  matrix[N_theta, K_theta] X_theta; // Design matrix for regressors on theta
  matrix[S, N_theta] G_theta; // Groups for theta
  matrix[N_gamma, K_gamma] X_gamma; // Design matrix for regressors on gamma
  matrix[S, N_gamma] G_gamma; // Groups for gamma
  real<lower=0> pa_2017; // Price of cattle in 2017

  // Prior hyperparams
  cov_matrix[K_theta] inv_Lambda_theta;
  vector[K_theta] mu_theta;
  real<lower=0> a_theta;
  real<lower=0> b_theta;

  cov_matrix[K_gamma] inv_Lambda_gamma;
  vector[K_gamma] mu_gamma;
  real<lower=0> a_gamma;
  real<lower=0> b_gamma;
}
generated quantities {
  // Priors
  real<lower=0> sigma_sq_theta = inv_gamma_rng(a_theta, b_theta);
  vector[K_theta] beta_theta = multi_normal_rng(mu_theta, sigma_sq_theta * inv_Lambda_theta);
  array[N_theta] real<lower=0> nabla_theta = lognormal_rng(X_theta * beta_theta, sqrt(sigma_sq_theta));

  real<lower=0> sigma_sq_gamma = inv_gamma_rng(a_gamma, b_gamma);
  vector[K_gamma] beta_gamma = multi_normal_rng(mu_gamma, sigma_sq_gamma * inv_Lambda_gamma);
  array[N_gamma] real<lower=0> nabla_gamma = lognormal_rng(X_gamma * beta_gamma, sqrt(sigma_sq_gamma));

  // Grouped average
  vector<lower=0>[S] theta = (G_theta * to_vector(nabla_theta)) / pa_2017;
  vector<lower=0>[S] gamma = G_gamma * to_vector(nabla_gamma);
}
