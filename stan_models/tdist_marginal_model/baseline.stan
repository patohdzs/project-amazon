data {
  int<lower=0> S; // Number of sites
  int<lower=0> N_theta;
  int<lower=0> N_gamma;
  int<lower=0> K_theta; // Number of coefficients on theta
  int<lower=0> K_gamma; // Number of coefficients on gamma

  matrix[N_theta, K_theta] X_theta; // Design matrix for regressors on theta
  matrix[S, N_theta] G_theta; // Groups for theta

  matrix[N_gamma, K_gamma] X_gamma; // Design matrix for regressors on gamma
  matrix[S, N_gamma] G_gamma; // Groups for gamma
  real<lower=0> pa_2017; // Price of cattle in 2017

  // Prior hyperparams
  cov_matrix[K_theta] inv_Q_theta;
  vector[K_theta] m_theta;
  real<lower=0> a_theta;
  real<lower=0> b_theta;

  cov_matrix[K_gamma] inv_Q_gamma;
  vector[K_gamma] m_gamma;
  real<lower=0> a_gamma;
  real<lower=0> b_gamma;
}
transformed data {
  cov_matrix[K_theta] Sigma_theta = (b_theta / a_theta) * inv_Q_theta;
  cov_matrix[K_gamma] Sigma_gamma = (b_gamma / a_gamma) * inv_Q_gamma;
}
generated quantities {
  // Priors
  beta_theta ~ multi_student_t_rng(2 * a_theta, m_theta, Sigma_theta);
  beta_gamma ~ multi_student_t_rng(2 * a_gamma, m_gamma, Sigma_gamma);

  // Grouped average
  vector<lower=0>[S] theta = (G_theta * exp(X_theta * beta_theta)) / pa_2017;
  vector<lower=0>[S] gamma = G_gamma * exp(X_gamma * beta_gamma);
}
