functions {
  real log_value(vector gamma, vector theta, int T, int S,
                            real alpha, matrix sol_val_X, vector sol_val_Ua,
                            matrix sol_val_Up, vector zbar_2017,
                            vector forest_area_2017, vector alpha_p_Adym,
                            matrix Bdym, vector ds_vect, real zeta, real xi,
                            real kappa, real pa, real pf) {
    // Carbon captured at time zero
    real x0 = (gamma' * forest_area_2017);

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

    matrix[S, T ] z_shifted_X;
    vector[S] scl = pa * theta - pf * kappa;
    for (j in 1 : (T )) {
      z_shifted_X[ : , j] = sol_val_X[1 : S, j+1] .* scl;
    }

    // Adjustment costs
    real term_1 = -sum(ds_vect[1 : T] .* sol_val_Ua) * zeta / 2.0;

    // Value of emissions absorbed
    real term_2 = sum(ds_vect[1 : T] .* (X_dym[2 : (T + 1)] - X_dym[1 : T]))
                  * pf;

    // Value of cattle output minus cost of emissions
    real term_3 = sum(ds_vect[1 : T] .* (rep_row_vector(1.0, S) * z_shifted_X)');

    // Overall objective value
    real obj_val = term_1 + term_2 + term_3;
    real log_density_val = -1.0 / xi * obj_val;

    return log_density_val;
  }
}
data {
  int<lower=0> T; // Time horizon
  int<lower=0> S; // Number of sites
  real<lower=0> alpha; // Mean reversion coefficient
  matrix[S + 2, T + 1] sol_val_X; // Sate trajectories
  vector[T] sol_val_Ua; // Squared total control adjustments; dimensions T x 1
  matrix[S, T] sol_val_Up; // U control
  vector[S] zbar_2017; // z_bar in 2017
  vector[S] forest_area_2017; // forrest area in 2017
  vector[T] alpha_p_Adym;
  matrix[T, T] Bdym;
  vector[T + 1] ds_vect; // Time discounting vector
  real zeta; // Penalty on adjustment costs
  real xi; // Penalty on prior-posterior KL div
  real kappa; // Effect of cattle farming on emissions
  real<lower=0> pa; // Price of cattle output
  real<lower=0> pf; // Price of carbon emission transfers




  // theta parameters

  int<lower=0> N_theta; // Number of observations on theta
  int<lower=0> K_theta; // Number of coefficients on theta
  int<lower=0> M_theta; // Number of large groups
  matrix[N_theta, K_theta] X_theta; // Design matrix for regressors on theta
  matrix[S, N_theta] SG_theta; // Site Groups for theta

  vector[K_theta] theta_mean;
  cov_matrix[K_theta] theta_vcov;
  vector[M_theta] V_theta_mean;
  vector[M_theta] V_theta_var;
  array[N_theta] int G_theta;

  real<lower=0> pa_2017; // Price of cattle in 2017



  // gamma parameters

  int<lower=0> N_gamma; // Number of observations on gamma
  int<lower=0> K_gamma; // Number of coefficients on gamma
  int<lower=0> M_gamma; // Number of large groups
  matrix[N_gamma, K_gamma] X_gamma; // Design matrix for regressors on gamma

  vector[K_gamma] gamma_mean;
  cov_matrix[K_gamma] gamma_vcov;
  vector[M_gamma] V_gamma_mean;
  vector[M_gamma] V_gamma_var;
  array[N_gamma] int G_gamma;


}


transformed data {
  matrix[K_theta, K_theta] L_theta = cholesky_decompose(theta_vcov);
  matrix[K_gamma, K_gamma] L_gamma = cholesky_decompose(gamma_vcov);  
}


parameters {
  vector[K_theta] alpha_theta;
  vector[M_theta] Vj_theta;

  vector[K_gamma] alpha_gamma;
  vector[M_gamma] Vj_gamma;

}
transformed parameters {
  // Coefs
  vector[K_theta] beta_theta = theta_mean +  L_theta * alpha_theta;

  vector[K_gamma] beta_gamma = gamma_mean +  L_gamma * alpha_gamma;

  // Grouped average
  vector<lower=0>[S] theta = (SG_theta * exp(X_theta * beta_theta + Vj_theta[G_theta])) / pa_2017;
  
  vector<lower=0>[S] gamma = exp(X_gamma * beta_gamma + Vj_gamma[G_gamma]);

}
model {
  // Hierarchical priors

  alpha_theta ~ std_normal();

  for (j in 1:M_theta) {
  Vj_theta[j] ~ normal(V_theta_mean[j],sqrt(V_theta_var[j]));
  }


  alpha_gamma ~ std_normal();

  for (j in 1:M_gamma) {
  Vj_gamma[j] ~ normal(V_gamma_mean[j],sqrt(V_gamma_var[j]));
  }



  // Value function
  target += log_value(gamma, theta, T, S, alpha, sol_val_X,
                                 sol_val_Ua, sol_val_Up, zbar_2017,
                                 forest_area_2017, alpha_p_Adym, Bdym,
                                 ds_vect, zeta, xi, kappa, pa, pf);
}


  
