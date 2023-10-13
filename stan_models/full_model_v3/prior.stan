data {
  int<lower=0> S; // Number of sites
  vector[N_theta] y_theta;
  matrix[N_theta, K_theta] X_theta; // Design matrix for regressors on theta
  matrix[S, N_theta] G_theta; // Groups for theta

  vector[N_gamma] y_gamma;
  matrix[N_gamma, K_gamma] X_gamma; // Design matrix for regressors on gamma
  matrix[S, N_gamma] G_gamma; // Groups for gamma
  real<lower=0> pa_2017; // Price of cattle in 2017
}
transformed data {
  cov_matrix[K_theta] inv_Lambda_theta = inverse_spd(X_theta' * X_theta);
  vector[K_theta] mu_theta = inv_Lambda_theta * X_theta' y_theta;
  real<lower=0> a_theta = N_theta / 2;
  real<lower=0> b_theta = (y_theta' * y_theta - mu_theta' * X' * X * mu_theta) / 2;

  cov_matrix[K_gamma] inv_Lambda_gamma = inverse_spd(X_gamma' * X_gamma);
  vector[K_gamma] mu_gamma = inv_Lambda_gamma * X_gamma' * y_gamma ;
  real<lower=0> a_gamma = N_gamma / 2;
  real<lower=0> b_gamma = (y_gamma' * y_gamma - mu_gamma' * X_gamma' * X_gamma * mu_gamma) / 2;
}
generated quantities {
  // Priors
  real<lower=0> sigma_sq_theta = inv_gamma_rng(a_theta, b_theta);
  vector[K_theta] beta_theta = multi_normal_rng(mu_theta, sigma_sq_theta * inv_Lambda_theta);

  real<lower=0> sigma_sq_gamma = inv_gamma_rng(a_gamma, b_gamma);
  vector[K_gamma] beta_gamma = multi_normal_rng(mu_gamma, sigma_sq_gamma * inv_Lambda_gamma);

  // Grouped average
  vector<lower=0>[S] theta = (G_theta * exp(X_theta * beta_theta)) / pa_2017;
  vector<lower=0>[S] gamma = G_gamma * exp(X_gamma * beta_gamma);
}
