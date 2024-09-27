data {
  int<lower=0> S; // Number of sites
  int<lower=0> N_theta;
  int<lower=0> K_theta; // Number of coefficients on theta


  matrix[N_theta, K_theta] X_theta; // Design matrix for regressors on theta
  matrix[S, N_theta] G_theta; // Groups for theta


  real<lower=0> pa_2017; // Price of cattle in 2017

  // Prior hyperparams
  cov_matrix[K_theta] inv_Q_theta;
  vector[K_theta] m_theta;
  real<lower=0> a_theta;
  real<lower=0> b_theta;

}
generated quantities {
  // Priors
  real<lower=0> sigma_sq_theta = inv_gamma_rng(a_theta, b_theta);
  vector[K_theta] beta_theta = multi_normal_rng(m_theta,
                                                sigma_sq_theta * inv_Q_theta);


  //Grouped average
  vector<lower=0>[S] theta = (G_theta * exp(X_theta * beta_theta)) / pa_2017;

}
