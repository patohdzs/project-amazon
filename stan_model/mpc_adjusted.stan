functions {
  real log_value(vector gamma, vector theta, int T, int S, int J, real alpha,
                array[] matrix Z, array[] matrix U, array[] matrix V, array[] matrix X,
                 vector zbar_2017,  vector alpha_p_Adym, matrix Bdym,
                 vector ds_vect, real zeta_u, real zeta_v, real xi,
                 real kappa, matrix pa, real pe, vector prob) {
    // Carbon captured at time zero

    real log_density_val = 0;


    for (j in 1:J) {
      
      matrix[S, T + 1] Z_j = Z[j];
      matrix[S, T] U_j = U[j];
      matrix[S, T] V_j = V[j];
      matrix[S, T + 1] X_j = X[j];
      vector[T+1] X_agg;
      for (t in 1:T+1) {
        X_agg[t] = sum(X_j[:,t]);
      }

      // Compute forest area in each time period for state j
      matrix[S, T] forest_area;
      for (t in 1:T) {
        forest_area[:, t] = zbar_2017 - Z_j[1:S, t];
      }


      // Compute aggregate X dot for state j
      vector[T] Xdot_agg = (X_agg[2:(T + 1)] - X_agg[1:T]);

      // Compute aggregate Z for state j
      vector[T] Z_agg = (rep_row_vector(1.0, S) * Z_j[:, 2:T + 1])';

      // Value of emissions absorbed for state j
      real term_1 = -pe * sum(ds_vect .* (kappa * Z_agg - Xdot_agg));

      // Value of agricultural output for state j
      vector[T + 1] agri_output = pa[ :, j ] .* (theta' * Z_j)';
      real term_2 = sum(ds_vect .* agri_output[2:T + 1]);

      // Value of adjustment costs for state j
      vector[T] U_agg = (rep_row_vector(1.0, S) * U_j)';
      vector[T] V_agg = (rep_row_vector(1.0, S) * V_j)';

      vector[T] w = zeta_u / 2.0 * (U_agg .* U_agg)
                    + zeta_v / 2.0 * (V_agg .* V_agg);

      real term_3 = -sum(ds_vect .* w);

      // Overall objective value for state j
      real obj_val_j = term_1 + term_2 + term_3;

      // Add weighted contribution of state j to the total log density value
      log_density_val += prob[j] * (-1.0 / xi * obj_val_j);
    }

    return log_density_val;
  }
}
data {
  // Planner problem
  int<lower=0> T; // Time horizon
  int<lower=0> S; // Number of sites
  int<lower=0> J; // Number of states (J=8)
  array[J] matrix[S, T + 1] Z; // Agricultural area state
  array[J] matrix[S, T + 1] X;
  array[J] matrix[S, T] U; // U control
  array[J] matrix[S, T] V; // V control
  vector[S] zbar_2017; // z_bar in 2017
  vector[T] alpha_p_Adym;
  matrix[T, T] Bdym;
  vector[T] ds_vect; // Time discounting vector
  real<lower=0> alpha; // Mean-reversion coefficient
  real zeta_u; // Penalty on adjustment costs
  real zeta_v; // Penalty on adjustment costs
  real xi; // Penalty on prior-posterior KL div
  real kappa; // Effect of cattle farming on emissions
  matrix[T+1, J] pa; // Price of cattle output
  real<lower=0> pe; // Price of carbon emission transfers
  real<lower=0> pa_current;
  vector[S] gamma;
  vector[S] theta;

}

parameters {
  simplex[2] p_ll;  // Transition probabilities for state LL
  simplex[2] p_hh;  // Transition probabilities for state HH
}

transformed parameters {
  vector[8] prob;  // Probabilities for the 8 states

  if (pa_current == 1) {
    prob[1] = p_ll[1]^3;
    prob[2] = p_ll[1]^2 * (1 - p_ll[1]);
    prob[3] = p_ll[1] * (1 - p_ll[1]) * (1 - p_hh[1]);
    prob[4] = p_ll[1] * (1 - p_ll[1]) * p_hh[1];
    prob[5] = (1 - p_ll[1]) * (1 - p_hh[1]) * p_ll[1];
    prob[6] = (1 - p_ll[1]) * (1 - p_hh[1]) * (1 - p_ll[1]);
    prob[7] = (1 - p_ll[1]) * p_hh[1] * (1 - p_hh[1]);
    prob[8] = (1 - p_ll[1]) * p_hh[1] * p_hh[1];
  } else if (pa_current == 2) {
    prob[1] = (1 - p_hh[1]) * p_ll[1]^2;
    prob[2] = (1 - p_hh[1]) * p_ll[1] * (1 - p_ll[1]);
    prob[3] = (1 - p_hh[1]) * (1 - p_ll[1]) * (1 - p_hh[1]);
    prob[4] = (1 - p_hh[1]) * (1 - p_ll[1]) * p_hh[1];
    prob[5] = p_hh[1] * (1 - p_hh[1]) * p_ll[1];
    prob[6] = p_hh[1] * (1 - p_hh[1]) * (1 - p_ll[1]);
    prob[7] = p_hh[1] * p_hh[1] * (1 - p_hh[1]);
    prob[8] = p_hh[1] * p_hh[1] * p_hh[1];
  }
}


model {
  p_ll[1] ~ beta(14.8, 6.2);  // Mean = 0.706, Var = 0.01
  p_hh[1] ~ beta(21.4, 4.4);  // Mean = 0.829, Var = 0.01
  // p_ll[1] ~ uniform(0.705, 0.707);  // Tight range around 0.706
  // p_hh[1] ~ uniform(0.828, 0.830);  // Tight range around 0.829

  // Value function
  target += log_value(gamma, theta, T, S, J,alpha, Z, U, V, X, zbar_2017,
                      alpha_p_Adym, Bdym, ds_vect, zeta_u,
                      zeta_v, xi, kappa, pa, pe, prob);
}
