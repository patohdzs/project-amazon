functions {
  real planner_obj(int T, int S, vector theta, vector gamma,
                    matrix z, matrix u, matrix v, vector w,
                    vector z_bar, vector forest_area_2017,
                    real dt, real pe, real pa,
                    real alpha, real delta, real kappa,
                    real zeta, real xi) {

    // Initial condition for site emissions
    matrix[T + 1, S] x;
    x[1,:] = (gamma .* forest_area_2017)';

    // Computing site emissions
    for (t in 1 : T) {
      x[t + 1,:] = x[t,:] + (-gamma' .* u[t, :] - alpha * (x[t,:] - gamma' .* (z_bar' - z[t, :]))) * dt;
    }

    // Computing result
    real result = 0;

    for (t in 1 : 50) {
      real discount_factor = exp(-delta * (t * dt - dt));
      real emissions = -pe * sum(kappa * z[t, :] - (x[t + 1, :] - x[t, :]) / dt);
      real output = pa * sum(theta .* z[t, :]');
      real adj_costs = -zeta / 2 * w[t]^2;

      result += discount_factor * (emissions + output + adj_costs) * dt;
    }

    return (-1.0 / xi) * result;
  }
}
data {
  int<lower=0> T;
  int<lower=0> S;
  matrix[T + 1, S] z;
  matrix[T + 1, S] u;
  matrix[T + 1, S] v;
  vector[T + 1] w;
  real delta;
  real<lower=0> dt;
  real<lower=0> pe;
  real<lower=0> pa;

  vector[S] forest_area_2017;
  vector[S] z_bar;

  real<lower=0> alpha;
  real kappa;
  real zeta;
  real xi;

  int<lower=0> N_theta;
  int<lower=0> N_gamma;
  int<lower=0> K_theta;
  int<lower=0> K_gamma;

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
  matrix[K_theta, K_theta] L_theta = cholesky_decompose(inv_Q_theta);
  matrix[K_gamma, K_gamma] L_gamma = cholesky_decompose(inv_Q_gamma);
}
parameters {
  real<lower=0> sigma_sq_theta; // Variance of log_theta
  vector[K_theta] alpha_theta;

  real<lower=0> sigma_sq_gamma; // Variance of log_gamma
  vector[K_gamma] alpha_gamma;
}
transformed parameters {
  // Coefs
  vector[K_theta] beta_theta = m_theta + sqrt(sigma_sq_theta) * L_theta * alpha_theta;
  vector[K_gamma] beta_gamma = m_gamma + sqrt(sigma_sq_gamma) * L_gamma * alpha_gamma;

  // Grouped average
  vector<lower=0>[S] theta = (G_theta * exp(X_theta * beta_theta + sigma_sq_theta/2)) / pa_2017;
  vector<lower=0>[S] gamma = G_gamma * exp(X_gamma * beta_gamma + sigma_sq_gamma/2);
}
model {
  // Hierarchical priors
  sigma_sq_theta ~ inv_gamma(a_theta, b_theta);
  sigma_sq_gamma ~ inv_gamma(a_gamma, b_gamma);

  alpha_theta ~ std_normal();
  alpha_gamma ~ std_normal();

  // Value function
  target += planner_obj(T, S, theta, gamma,
                     z, u, v, w, z_bar,
                     forest_area_2017, dt,
                     pe, pa,alpha, delta, kappa,
                     zeta, xi);
}
