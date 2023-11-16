functions {
  real log_value(vector gamma, vector theta, int T, int S,
                            real alpha, matrix sol_val_X, vector sol_val_Ua,
                            matrix sol_val_Up, vector zbar_2017,
                            vector forestArea_2017_ha, vector alpha_p_Adym,
                            matrix Bdym, vector ds_vect, real zeta, real xi,
                            real kappa, real pa, real pf) {
    // Carbon captured at time zero
    real x0 = (gamma' * forestArea_2017_ha);

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
  real planner_obj(int T, int S, vector theta,
                    matrix z, matrix x, vector w, real delta,
                    real dt, real pe, real pa,
                    real kappa, real zeta, real xi) {
    real result = 0;

    for (t in 1 : T) {
      if (t < T) {
        real exp_term = exp(-delta * (t * dt - dt));
        real pe_term = -pe * sum(kappa * z[t, :] - (x[t + 1, :] - x[t, :]) / dt);
        real pa_term = pa * sum(theta .* z[t, :]');
        real zeta_term = -zeta / 2 * w[t]^2;

        result += exp_term * (pe_term + pa_term + zeta_term) * dt;
      }
    }

    return (-1.0 / xi) * result;
  }
}
data {
  int<lower=0> T;
  int<lower=0> S;
  matrix[T + 1, S] z;
  matrix[T + 1, S] x;
  vector[T + 1] w;
  real delta;
  real<lower=0> dt;
  real<lower=0> pe;
  real<lower=0> pa;
  real kappa;
  real zeta;
  real xi;

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
transformed data {
  matrix[K_theta, K_theta] L_theta = cholesky_decompose(inv_Q_theta);
}
parameters {
  real<lower=0> sigma_sq_theta; // Variance of log_theta
  vector[K_theta] alpha_theta;
}
transformed parameters {
  // Coefs
  vector[K_theta] beta_theta = m_theta + sqrt(sigma_sq_theta) * L_theta * alpha_theta;

  // Grouped average
  vector<lower=0>[S] theta = (G_theta * exp(X_theta * beta_theta + sigma_sq_theta/2)) / pa_2017;
}
model {
  // Hierarchical priors
  sigma_sq_theta ~ inv_gamma(a_theta, b_theta);

  alpha_theta ~ std_normal();

  // Value function
  target += planner_obj(T, S, theta, z, x, w, delta,
                        dt, pe, pa, kappa, zeta, xi);
}
