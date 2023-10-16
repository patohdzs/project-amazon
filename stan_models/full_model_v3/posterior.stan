functions {
  real log_density_function(vector gamma, vector theta, int T, int S,
                            real alpha, matrix sol_val_X, vector sol_val_Ua,
                            matrix sol_val_Up, vector zbar_2017,
                            vector forestArea_2017_ha, real norm_fac,
                            vector alpha_p_Adym, matrix Bdym, vector ds_vect,
                            real zeta, real xi, real kappa, real pa, real pf) {
    // Carbon captured at time zero
    real x0 = (gamma' * forestArea_2017_ha) / norm_fac;

    // Aggregate carbon captured
    vector[T] X_zero = x0 * rep_vector(1.0, T);

    matrix[S, T] shifted_X;
    for (j in 1 : T) {
      shifted_X[ : , j] = zbar_2017 - sol_val_X[1 : S, j];
    }
    row_vector[T] omega = gamma' * ((alpha * shifted_X) - sol_val_Up);

    vector[T + 1] X_dym;
    X_dym[1] = x0;
    X_dym[2 : (T + 1)] = alpha_p_Adym .* X_zero + Bdym * omega';

    matrix[S, T + 1] z_shifted_X;
    vector[S] scl = pa * theta - pf * kappa;
    for (j in 1 : (T + 1)) {
      z_shifted_X[ : , j] = sol_val_X[1 : S, j] .* scl;
    }

    // Adjustment costs
    real term_1 = -sum(ds_vect[1 : T] .* sol_val_Ua) * zeta / 2.0;

    // Value of emissions absorbed
    real term_2 = sum(ds_vect[1 : T] .* (X_dym[2 : (T + 1)] - X_dym[1 : T]))
                  * pf;

    // Value of cattle output minus cost of emissions
    real term_3 = sum(ds_vect .* (rep_row_vector(1.0, S) * z_shifted_X)');

    // Overall objective value
    real obj_val = term_1 + term_2 + term_3;
    real log_density_val = -1.0 / xi * obj_val;

    return log_density_val;
  }
}
data {
  int<lower=0> T; // Time horizon
  int<lower=0> S; // Number of sites
  real<lower=0> norm_fac; // Normalization factor
  real<lower=0> alpha; // Mean reversion coefficient
  matrix[S + 2, T + 1] sol_val_X; // Sate trajectories
  vector[T] sol_val_Ua; // Squared total control adjustments; dimensions T x 1
  matrix[S, T] sol_val_Up; // U control
  vector[S] zbar_2017; // z_bar in 2017
  vector[S] forestArea_2017_ha; // forrest area in 2017
  vector[T] alpha_p_Adym;
  matrix[T, T] Bdym;
  vector[T + 1] ds_vect; // Time discounting vector
  real zeta; // Penalty on adjustment costs
  real xi; // Penalty on prior-posterior KL div
  real kappa; // Effect of cattle farming on emissions
  real<lower=0> pa; // Price of cattle output
  real<lower=0> pf; // Price of carbon emission transfers

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
  cov_matrix[K_theta] inv_Lambda_theta;
  vector[K_theta] mu_theta;
  real<lower=0> a_theta;
  real<lower=0> b_theta;

  cov_matrix[K_gamma] inv_Lambda_gamma;
  vector[K_gamma] mu_gamma;
  real<lower=0> a_gamma;
  real<lower=0> b_gamma;
}
parameters {
  real<lower=0> sigma_sq_theta; // Variance of log_theta
  vector[K_theta] beta_theta; // Coefficients on theta

  real<lower=0> sigma_sq_gamma; // Variance of log_gamma
  vector[K_gamma] beta_gamma; // Coefficients on gamma
}
transformed parameters {
  // Grouped average
  vector<lower=0>[S] theta = (G_theta * exp(X_theta * beta_theta)) / pa_2017;
  vector<lower=0>[S] gamma = G_gamma * exp(X_gamma * beta_gamma);
}
model {
  // Hierarchical priors
  sigma_sq_theta ~ inv_gamma(a_theta, b_theta);
  beta_theta ~ multi_normal(mu_theta, sigma_sq_theta * inv_Lambda_theta);

  sigma_sq_gamma ~ inv_gamma(a_gamma, b_gamma);
  beta_gamma ~ multi_normal(mu_gamma, sigma_sq_gamma * inv_Lambda_gamma);

  // Value function
  target += log_density_function(gamma, theta, T, S, alpha, sol_val_X,
                                 sol_val_Ua, sol_val_Up, zbar_2017,
                                 forestArea_2017_ha, norm_fac, alpha_p_Adym,
                                 Bdym, ds_vect, zeta, xi, kappa, pa, pf);
}
