import time

import numpy as np
from cmdstanpy import CmdStanModel
import jax.numpy as jnp
from jax import jit
import jax
from jax.scipy.linalg import cholesky
import numpyro
import numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS
from functools import partial

from pysrc.sampling import baseline

from ..optimization import gams, gurobi
from ..sampling import gamma_adj_reg_data, theta_adj_reg_data
from ..sampling.baseline import baseline_hyperparams
from ..services.data_service import load_site_data
from ..services.file_service import stan_model_path

@partial(jax.jit, static_argnums=(2,3))
def log_value_jax(gamma, theta, T, S, alpha, sol_val_X, sol_val_Ua, sol_val_Up, zbar_2017, forest_area_2017, alpha_p_Adym, Bdym, ds_vect, zeta, xi, kappa, pa, pf):
    # Carbon captured at time zero
    x0 = jnp.dot(gamma, forest_area_2017)

    # Aggregate carbon captured
    X_zero = x0 * jnp.ones(T)

    # Initialize shifted_X with zeros and then fill
    shifted_X = jnp.zeros((S, T))
    for j in range(T):  # Adjusted for zero-based indexing
        shifted_X = shifted_X.at[:, j].set(zbar_2017 - sol_val_X[:S, j])

    omega = jnp.dot(gamma, (alpha * shifted_X) - sol_val_Up)

    # Initialize X_dym and fill
    X_dym = jnp.zeros(T + 1)
    X_dym = X_dym.at[0].set(x0)
    X_dym = X_dym.at[1:].set(alpha_p_Adym * X_zero + jnp.dot(Bdym, omega.T))

    # Initialize z_shifted_X with zeros and then fill
    z_shifted_X = jnp.zeros((S, T))
    scl = pa * theta - pf * kappa
    for j in range(T):  # Adjusted for zero-based indexing
        z_shifted_X = z_shifted_X.at[:, j].set(sol_val_X[:S, j] * scl)  # Adjusted for sol_val_X indexing

    # Adjustment costs
    term_1 = -jnp.sum(ds_vect[:T] * sol_val_Ua) * zeta / 2.0

    # Value of emissions absorbed
    term_2 = jnp.sum(ds_vect[:T] * (X_dym[1:] - X_dym[:-1])) * pf

    # Value of cattle output minus cost of emissions
    term_3 = jnp.sum(ds_vect[:T] * jnp.sum(z_shifted_X, axis=0))

    # Overall objective value
    obj_val = term_1 + term_2 + term_3
    log_density_val = -1.0 / xi * obj_val

    return log_density_val


def numpyro_model(T, S, alpha, sol_val_X, sol_val_Ua, sol_val_Up,sol_val_Um,sol_val_Z, zbar_2017, forest_area_2017, 
                  alpha_p_Adym, Bdym, ds_vect, zeta, xi, kappa, pa, pf, N_theta, N_gamma, 
                  K_theta, K_gamma, X_theta, G_theta, X_gamma, G_gamma, pa_2017, 
                  inv_Q_theta, m_theta, a_theta, b_theta, inv_Q_gamma, m_gamma, a_gamma, b_gamma):
    # Transformed data
    L_theta = cholesky(inv_Q_theta)
    L_gamma = cholesky(inv_Q_gamma)

    # Parameters
    sigma_sq_theta = numpyro.sample('sigma_sq_theta', dist.InverseGamma(a_theta, b_theta))
    alpha_theta = numpyro.sample('alpha_theta', dist.Normal(0, 1).expand([K_theta]))
    
    sigma_sq_gamma = numpyro.sample('sigma_sq_gamma', dist.InverseGamma(a_gamma, b_gamma))
    alpha_gamma = numpyro.sample('alpha_gamma', dist.Normal(0, 1).expand([K_gamma]))
    
    # Transformed parameters
    beta_theta = m_theta + jnp.sqrt(sigma_sq_theta) * jnp.dot(L_theta, alpha_theta)
    beta_gamma = m_gamma + jnp.sqrt(sigma_sq_gamma) * jnp.dot(L_gamma, alpha_gamma)
    beta_theta= numpyro.deterministic("beta_theta", beta_theta)
    beta_gamma= numpyro.deterministic("beta_gamma", beta_gamma)
    
    theta = jnp.exp(jnp.dot(X_theta, beta_theta)) / pa_2017
    gamma = jnp.exp(jnp.dot(X_gamma, beta_gamma))

    # Ensure non-negativity
    theta = numpyro.deterministic('theta', jnp.maximum(0, jnp.dot(G_theta, theta)))
    gamma = numpyro.deterministic('gamma', jnp.maximum(0, jnp.dot(G_gamma, gamma)))

    # Value function (log-probability contribution)
    log_val = log_value_jax(gamma, theta, T, S, alpha, sol_val_X, sol_val_Ua, sol_val_Up, 
                            zbar_2017, forest_area_2017, alpha_p_Adym, Bdym, ds_vect, 
                            zeta, xi, kappa, pa, pf)
    numpyro.factor('log_value', log_val)





def sample(
    model_name,
    xi,
    pe,
    pa,
    weight,
    num_sites,
    # Model parameters
    T,
    dt=1,
    alpha=0.045007414,
    delta=0.02,
    kappa=2.094215255,
    zeta=1.66e-4 * 1e9,  # use the same normalization factor
    pa_2017=44.9736197781184,
    # Optimizer
    optimizer="gurobi",
    # Sampling params
    max_iter=20000,
    tol=0.001,
    final_sample_size=5_000,
    **stan_kwargs,
):
    # Instantiate stan sampler
    sampler = CmdStanModel(
        stan_file=stan_model_path(model_name) / "adjusted.stan",
        cpp_options={"STAN_THREADS": "true"},
        force_compile=True,
    )

    # Load sites' data
    (
        zbar_2017,
        _,
        z_2017,
        forest_area_2017,
        _,
        site_theta_df,
        site_gamma_df,
        municipal_theta_df,
        municipal_gamma_df,
    ) = load_site_data(num_sites)

    # Set initial theta & gamma using baseline mean
    fit = baseline.sample(
        model_name=model_name, num_samples=final_sample_size, num_sites=num_sites
    )
    theta_vals = fit.stan_variable("theta").mean(axis=0)
    gamma_vals = fit.stan_variable("gamma").mean(axis=0)

    # Save starting params
    uncertain_vals = np.concatenate((theta_vals, gamma_vals)).copy()
    uncertain_vals_old = np.concatenate((theta_vals, gamma_vals)).copy()

    # Collected Ensembles over all iterations; dictionary indexed by iteration number
    collected_ensembles = {}
    coe_ensembles = {}

    # Track error over iterations
    uncertain_vals_tracker = [uncertain_vals_old.copy()]
    abs_error_tracker = []
    pct_error_tracker = []
    solution_tracker = []
    sampling_time_tracker = []
    fit_tracker = []

    # Results dictionary
    results = dict(
        num_sites=num_sites,
        tol=tol,
        T=T,
        dt=dt,
        delta_t=delta,
        alpha=alpha,
        kappa=kappa,
        pf=pe,
        pa=pa,
        xi=xi,
        zeta=zeta,
        final_sample_size=final_sample_size,
        weight=weight,
    )

    # Initialize error & iteration counter
    abs_error = np.infty
    pct_error = np.infty
    cntr = 0

    # Loop until convergence
    while cntr < max_iter and pct_error > tol:
        print(f"Optimization Iteration[{cntr+1}/{max_iter}]\n")

        # Flatten uncertain values
        uncertain_vals = np.asarray(uncertain_vals).flatten()

        # Unpacking uncertain values
        theta_vals = uncertain_vals[:num_sites].copy()
        gamma_vals = uncertain_vals[num_sites:].copy()

        print(f"Theta: {theta_vals}\n")
        print(f"Gamma: {gamma_vals}\n")

        # Computing carbon absorbed in start period
        x0_vals = gamma_vals * forest_area_2017

        # Choose optimizer
        if optimizer == "gurobi":
            solve_planner_problem = gurobi.solve_planner_problem

        elif optimizer == "gams":
            solve_planner_problem = gams.solve_planner_problem

        else:
            raise ValueError("Optimizer must be one of ['gurobi', 'gams']")

        trajectories = solve_planner_problem(
            T=T,
            theta=theta_vals,
            gamma=gamma_vals,
            x0=x0_vals,
            z0=z_2017,
            zbar=zbar_2017,
            dt=dt,
            pe=pe,
            pa=pa,
            alpha=alpha,
            delta=delta,
            kappa=kappa,
            zeta=zeta,
        )

        # Update trackers
        solution_tracker.append(trajectories)

        # HMC sampling
        print("Starting HMC sampling...\n")
        
        
        
        
        nuts_kernel = NUTS(numpyro_model)
        print("test2")
# Setup the MCMC run
        mcmc = MCMC(nuts_kernel, num_warmup=500, num_samples=4000)
        print("test3")
        # Run the MCMC
        mcmc.run(jax.random.PRNGKey(0),
                T=T, S=num_sites, 
                alpha=alpha, 
                zbar_2017=zbar_2017, 
                forest_area_2017=forest_area_2017, 
                zeta=zeta, 
                xi=xi, 
                kappa=kappa, 
                pa=pa, 
                pa_2017=pa_2017,
                pf=pe, 
                **trajectories,
                **_dynamics_matrices(T, dt, alpha, delta),
                **theta_adj_reg_data(num_sites, site_theta_df),
                **gamma_adj_reg_data(num_sites, site_gamma_df),
                **baseline_hyperparams(municipal_theta_df, "theta"),
                **baseline_hyperparams(municipal_gamma_df, "gamma"),
            )

# Get the samples
        samples = mcmc.get_samples()
        
        
        
        # import pickle

        # # Assuming `samples` is the dictionary of samples returned by mcmc.get_samples()
        # with open('mcmc_samples.pcl', 'wb') as file:
        #     pickle.dump(samples, file)

        # print("Samples saved to 'numpyro.pcl'.")

        theta_adj_samples = samples['theta']
        gamma_adj_samples = samples['gamma']
        theta_coe_adj_samples = samples["beta_theta"]
        gamma_coe_adj_samples = samples["beta_gamma"]
        print("theta",theta_adj_samples.shape)
        print("theta_coe",theta_coe_adj_samples.shape)
        # print("end here")
        # import sys; sys.exit()
        
        
        
        # model_data = dict(
        #     T=T,
        #     S=num_sites,
        #     alpha=alpha,
        #     zbar_2017=zbar_2017,
        #     forest_area_2017=forest_area_2017,
        #     zeta=zeta,
        #     xi=xi,
        #     kappa=kappa,
        #     pa=pa,
        #     pa_2017=pa_2017,
        #     pf=pe,
        #     **trajectories,
        #     **_dynamics_matrices(T, dt, alpha, delta),
        #     **theta_adj_reg_data(num_sites, site_theta_df),
        #     **gamma_adj_reg_data(num_sites, site_gamma_df),
        #     **baseline_hyperparams(municipal_theta_df, "theta"),
        #     **baseline_hyperparams(municipal_gamma_df, "gamma"),
        # )

        # # Sampling from adjusted distribution
        # sampling_time = time.time()
        # fit = sampler.sample(
        #     data=model_data,
        #     **stan_kwargs,
        # )
        # sampling_time = time.time() - sampling_time
        # print(f"Finished sampling! Elapsed Time: {sampling_time} seconds\n")
        # print(fit.diagnose())

        # # Update fit and sampling time trackers
        # fit_tracker.append(fit.summary())
        # sampling_time_tracker.append(sampling_time)

        # # Extract samples
        # theta_adj_samples = fit.stan_variable("theta")
        # gamma_adj_samples = fit.stan_variable("gamma")
        # theta_coe_adj_samples = fit.stan_variable("beta_theta")
        # gamma_coe_adj_samples = fit.stan_variable("beta_gamma")
        # print("theta",theta_adj_samples.shape)
        # print("end here")
        # import sys; sys.exit()
        
        uncertainty_adj_samples = np.concatenate(
            (theta_adj_samples, gamma_adj_samples), axis=1
        )

        uncertainty_coe_adj_samples = np.concatenate(
            (theta_coe_adj_samples, gamma_coe_adj_samples), axis=1
        )

        # Update ensemble/tracker
        collected_ensembles.update({cntr: uncertainty_adj_samples.copy()})
        coe_ensembles.update({cntr: uncertainty_coe_adj_samples.copy()})

        print(f"Parameters from last iteration: {uncertain_vals_old}\n")
        print(
            f"""Parameters from current iteration:
            {np.mean(uncertainty_adj_samples, axis=0)}\n"""
        )

        # Compute exponentially-smoothened new params
        uncertain_vals = (
            weight * np.mean(uncertainty_adj_samples, axis=0)
            + (1 - weight) * uncertain_vals_old
        )

        uncertain_vals_tracker.append(uncertain_vals.copy())
        print(f"Updated uncertain values: {uncertain_vals}\n")

        # Evaluate error for convergence check
        # The percentage difference are changed to absolute difference
        abs_error = np.max(np.abs(uncertain_vals_old - uncertain_vals))
        pct_error = np.max(
            np.abs(uncertain_vals_old - uncertain_vals) / uncertain_vals_old
        )

        abs_error_tracker.append(abs_error)
        pct_error_tracker.append(pct_error)

        print(
            f"""
            Iteration [{cntr+1:4d}]: Absolute Error = {abs_error},
            Percentage Error = {pct_error}
            """
        )

        # Exchange parameter values (for future weighting/update & error evaluation)
        uncertain_vals_old = uncertain_vals

        # Increase the counter
        cntr += 1

        # Update results directory
        results.update(
            {
                "cntr": cntr,
                "abs_error_tracker": np.asarray(abs_error_tracker),
                "pct_error_tracker": np.asarray(pct_error_tracker),
                "uncertain_vals_tracker": np.asarray(uncertain_vals_tracker),
                "sampling_time_tracker": sampling_time_tracker,
                "collected_ensembles": collected_ensembles,
                "solution_tracker": solution_tracker,
                "coe_ensembles": coe_ensembles,
                "fit_tracker": fit_tracker,
            }
        )

    # Sample (densly) the final distribution
    print("Terminated. Sampling the final distribution...\n")
    
    mcmc = MCMC(nuts_kernel, num_warmup=1000, num_samples=5000)

    # Run the MCMC
    mcmc.run(jax.random.PRNGKey(0),
            T=T, S=num_sites, 
            alpha=alpha, 
            zbar_2017=zbar_2017, 
            forest_area_2017=forest_area_2017, 
            zeta=zeta, 
            xi=xi, 
            kappa=kappa, 
            pa=pa, 
            pa_2017=pa_2017,
            pf=pe, 
            **trajectories,
            **_dynamics_matrices(T, dt, alpha, delta),
            **theta_adj_reg_data(num_sites, site_theta_df),
            **gamma_adj_reg_data(num_sites, site_gamma_df),
            **baseline_hyperparams(municipal_theta_df, "theta"),
            **baseline_hyperparams(municipal_gamma_df, "gamma"),
        )

# Get the samples
    final_samples = mcmc.get_samples()

    theta_adj_samples = final_samples['theta']
    gamma_adj_samples = final_samples['gamma']
    theta_coe_adj_samples = final_samples['beta_theta']
    gamma_coe_adj_samples = final_samples['beta_gamma']
    # stan_kwargs["iter_sampling"] = final_sample_size
    # fit = sampler.sample(
    #     data=model_data,
    #     **stan_kwargs,
    # )

    # Extract samples
    # theta_adj_samples = fit.stan_variable("theta")
    # gamma_adj_samples = fit.stan_variable("gamma")
    # theta_coe_adj_samples = fit.stan_variable("beta_theta")
    # gamma_coe_adj_samples = fit.stan_variable("beta_gamma")

    final_samples = np.concatenate((theta_adj_samples, gamma_adj_samples), axis=1)
    final_samples_coe = np.concatenate(
        (theta_coe_adj_samples, gamma_coe_adj_samples), axis=1
    )

    results.update({"final_sample": final_samples})
    results.update({"final_sample_coe": final_samples_coe})

    return results


def _dynamics_matrices(T, dt, alpha, delta):
    # Create dynamics matrices
    arr = np.cumsum(
        np.triu(np.ones((T, T))),
        axis=1,
    ).T
    Bdym = (1 - alpha) ** (arr - 1)
    Bdym[Bdym > 1] = 0.0
    Adym = np.arange(1, T + 1)
    alpha_p_Adym = np.power(1 - alpha, Adym)

    # Other placeholders!
    ds_vect = np.exp(-delta * np.arange(T + 1) * dt)
    ds_vect = np.reshape(ds_vect, (ds_vect.size, 1)).flatten()
    return {"alpha_p_Adym": alpha_p_Adym, "Bdym": Bdym, "ds_vect": ds_vect}
