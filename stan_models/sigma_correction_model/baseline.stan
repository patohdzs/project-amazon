data {
  int<lower=0> S; // Number of sites

  int<lower=0> N_theta;
  int<lower=0> N_gamma;

  int<lower=0> K_theta; // Number of coefficients on theta
  int<lower=0> K_gamma; // Number of coefficients on gamma

  matrix[N_theta, K_theta] X_theta; // Design matrix for regressors on theta
  matrix[N_gamma, K_gamma] X_gamma; // Design matrix for regressors on gamma

  matrix[S, N_theta] G_theta; // Groups for theta
  matrix[S, N_gamma] G_gamma; // Groups for gamma

  real<lower=0> pa_2017; // Price of cattle in 2017

  // Prior hyperparams
  vector[N_theta] w_theta;
  cov_matrix[K_theta] inv_Q_theta;
  vector[K_theta] m_theta;
  real<lower=0> a_theta;
  real<lower=0> b_theta;

  cov_matrix[K_gamma] inv_Q_gamma;
  vector[K_gamma] m_gamma;
  real<lower=0> a_gamma;
  real<lower=0> b_gamma;
}
generated quantities {
  // Priors
  real<lower=0> sigma_sq_theta = inv_gamma_rng(a_theta, b_theta);
  real<lower=0> sigma_sq_gamma = inv_gamma_rng(a_gamma, b_gamma);

  vector[K_theta] beta_theta = multi_normal_rng(m_theta, sigma_sq_theta * inv_Q_theta);
  vector[K_gamma] beta_gamma = multi_normal_rng(m_gamma, sigma_sq_gamma * inv_Q_gamma);

  // Preds
  vector<lower=0>[N_theta] eta_theta = exp(X_theta * beta_theta + (sigma_sq_theta / 2));
  vector<lower=0>[N_gamma] eta_gamma = exp(X_gamma * beta_gamma + (sigma_sq_gamma / 2));

  // Grouped average
  vector<lower=0>[S] theta = (G_theta * eta_theta) / pa_2017;
  vector<lower=0>[S] gamma = G_gamma * eta_gamma;
}
