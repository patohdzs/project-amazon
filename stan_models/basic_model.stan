functions {
  real log_density_function(vector gamma, vector theta, int T, int S,
                            real alpha, matrix sol_val_X, vector sol_val_Ua,
                            matrix sol_val_Up, vector zbar_2017,
                            vector forestArea_2017_ha, real norm_fac,
                            vector alpha_p_Adym, matrix Bdym, vector ds_vect,
                            real zeta, real xi, real kappa, real pa, real pf) {
    // Carbon captured at time zero
    vector[S] x0_vals;
    for (i in 1 : S) {
      x0_vals[i] = (gamma[i] * forestArea_2017_ha[i]) / norm_fac;
    }

    // Aggregate carbon captured
    vector[T] X_zero = sum(x0_vals) * to_vector(rep_vector(1.0, T));

    matrix[S, T] shifted_X;
    for (j in 1 : T) {
      shifted_X[ : , j] = zbar_2017 - sol_val_X[1 : S, j];
    }
    row_vector[T] omega = gamma' * ((alpha * shifted_X) - sol_val_Up);

    vector[T + 1] X_dym;
    X_dym[1] = sum(x0_vals);
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
  int<lower=0> T; // time horizon
  int<lower=0> S; // number of sites
  real<lower=0> norm_fac; // Normalization factor
  real<lower=0> alpha; //
  matrix[S + 2, T + 1] sol_val_X; // Dim
  vector[T] sol_val_Ua; // Squared total control adjustments; dimensions T x 1
  matrix[S, T] sol_val_Up; // U control; dimensions I x T
  vector[S] zbar_2017;
  vector[S] forestArea_2017_ha;
  vector[T] alpha_p_Adym;
  matrix[T, T] Bdym;
  vector[T + 1] ds_vect;
  real zeta;
  real xi;
  real kappa;
  real<lower=0> pa;
  real<lower=0> pf;
}
parameters {
  vector[S] gamma; // gamma
  vector[S] theta; // theta
}
model {
  target += log_density_function(gamma, theta, T, S, alpha, sol_val_X,
                                 sol_val_Ua, sol_val_Up, zbar_2017,
                                 forestArea_2017_ha, norm_fac, alpha_p_Adym,
                                 Bdym, ds_vect, zeta, xi, kappa, pa, pf);
}
