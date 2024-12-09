functions {
  real log_value(vector gamma, vector theta, int T, int S, real alpha,
                 matrix Z, matrix U, matrix V, vector zbar_2017,
                 vector forest_area_2017, vector alpha_p_Adym, matrix Bdym,
                 vector ds_vect, real zeta_u, real zeta_v, real xi,
                 real kappa, real pa, real pe) {
    // Carbon captured at time zero
    real x0 = (gamma' * forest_area_2017);

    // Compute forest area in each time period
    matrix[S, T] forest_area;
    for (t in 1 : T) {
      forest_area[ : , t] = zbar_2017 - Z[1 : S, t];
    }

    // Aggregate carbon captured
    vector[T] X_zero = x0 * rep_vector(1.0, T);

    // Compute stock of carbon (X)
    row_vector[T] omega = gamma' * ((alpha * forest_area) - U);
    vector[T + 1] X;
    X[1] = x0;
    X[2 : (T + 1)] = alpha_p_Adym .* X_zero + Bdym * omega';

    // Compute aggregate X dot
    vector[T] Xdot_agg = (X[2 : (T + 1)] - X[1 : T]);

    // Compute aggregate Z
    vector[T] Z_agg = (rep_row_vector(1.0, S) * Z[ : , 2 : T + 1])';

    // Value of emissions absorbed
    real term_1 = -pe * sum(ds_vect .* (kappa * Z_agg - Xdot_agg));

    // Value of agricultural output
    vector[T + 1] agri_output = pa * (theta' * Z)';
    real term_2 = sum(ds_vect .* agri_output[2 : T + 1]);

    // Value of adjustment costs
    vector[T] U_agg = (rep_row_vector(1.0, S) * U)';
    vector[T] V_agg = (rep_row_vector(1.0, S) * V)';

    vector[T] w = zeta_u / 2.0 * (U_agg .* U_agg)
                  + zeta_v / 2.0 * (V_agg .* V_agg);

    real term_3 = -sum(ds_vect .* w);

    // Overall objective value
    real obj_val = term_1 + term_2 + term_3;
    real log_density_val = -1.0 / xi * obj_val;

    return log_density_val;
  }
}
data {
  // Planner problem
  int<lower=0> T; // Time horizon
  int<lower=0> S; // Number of sites
  matrix[S, T + 1] Z; // Agricultural area state
  matrix[S, T] U; // U control
  matrix[S, T] V; // V control
  vector[S] zbar_2017; // z_bar in 2017
  vector[S] forest_area_2017; // forest area in 2017
  vector[T] alpha_p_Adym;
  matrix[T, T] Bdym;
  vector[T] ds_vect; // Time discounting vector
  real<lower=0> alpha; // Mean-reversion coefficient
  real zeta_u; // Penalty on adjustment costs
  real zeta_v; // Penalty on adjustment costs
  real xi; // Penalty on prior-posterior KL div
  real kappa; // Effect of cattle farming on emissions
  real<lower=0> pa; // Price of cattle output
  real<lower=0> pe; // Price of carbon emission transfers

  // Theta parameters
  int<lower=0> N_theta; // Number of observations on theta
  int<lower=0> K_theta; // Number of coefficients on theta
  int<lower=0> M_theta; // Number of large groups
  matrix[N_theta, K_theta] X_theta; // Design matrix for regressors on theta
  matrix[S, N_theta] SG_theta; // Site Groups for theta

  vector[K_theta] beta_theta_mean; // Baseline distribution
  cov_matrix[K_theta] beta_theta_vcov;
  vector[M_theta] V_theta_mean;
  vector[M_theta] V_theta_var;
  array[N_theta] int G_theta; // Large group indicator

  real<lower=0> pa_2017; // Price of cattle in 2017

  // Gamma parameters
  int<lower=0> N_gamma; // Number of observations on gamma
  int<lower=0> K_gamma; // Number of coefficients on gamma
  int<lower=0> M_gamma; // Number of large groups
  matrix[N_gamma, K_gamma] X_gamma; // Design matrix for regressors on gamma

  vector[K_gamma] beta_gamma_mean; // Baseline distribution
  cov_matrix[K_gamma] beta_gamma_vcov;
  vector[M_gamma] V_gamma_mean;
  vector[M_gamma] V_gamma_var;
  array[N_gamma] int G_gamma; // Large group indicator
}
transformed data {
  matrix[K_theta, K_theta] L_theta = cholesky_decompose(beta_theta_vcov);
  matrix[K_gamma, K_gamma] L_gamma = cholesky_decompose(beta_gamma_vcov);
}
parameters {
  vector[K_theta] alpha_theta;
  vector[M_theta] Vj_theta;

  vector[K_gamma] alpha_gamma;
  vector[M_gamma] Vj_gamma;
}
transformed parameters {
  // Coefs
  vector[K_theta] beta_theta = beta_theta_mean + L_theta * alpha_theta;

  vector[K_gamma] beta_gamma = beta_gamma_mean + L_gamma * alpha_gamma;

  // Grouped average
  vector<lower=0>[S] theta = (SG_theta
                              * exp(X_theta * beta_theta + Vj_theta[G_theta]))
                             / pa_2017;

  vector<lower=0>[S] gamma = exp(X_gamma * beta_gamma + Vj_gamma[G_gamma]);
}
model {
  // Hierarchical priors

  alpha_theta ~ std_normal();

  for (j in 1 : M_theta) {
    Vj_theta[j] ~ normal(V_theta_mean[j], sqrt(V_theta_var[j]));
  }

  alpha_gamma ~ std_normal();

  for (j in 1 : M_gamma) {
    Vj_gamma[j] ~ normal(V_gamma_mean[j], sqrt(V_gamma_var[j]));
  }

  // Value function
  target += log_value(gamma, theta, T, S, alpha, Z, U, V, zbar_2017,
                      forest_area_2017, alpha_p_Adym, Bdym, ds_vect, zeta_u,
                      zeta_v, xi, kappa, pa, pe);
}
