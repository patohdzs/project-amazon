import os
import pickle
import time

import numpy as np
import stan
from optimization.casadi import solve_outer_optimization_problem
from sampling import gamma_reg_data, theta_reg_data
from services.data_service import load_site_data
from services.file_service import stan_model_path


def sample(
    model_name,
    output_dir,
    xi,
    pf,
    pa,
    weight,
    site_num,
    T,
    N=200,
    norm_fac=1e9,
    delta_t=0.02,
    alpha=0.045007414,
    kappa=2.094215255,
    zeta=1.66e-4 * 1e9,  # use the same normalization factor
    pa_2017=44.9736197781184,
    # Sampling params
    max_iter=20000,
    tol=0.001,
    sample_size=1000,
    final_sample_size=5_000,
    num_chains=8,
    num_warmup=500,
):
    # Create the output directory
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # Read model code
    with open(stan_model_path(model_name) / "posterior.stan") as f:
        model_code = f.read()

    # Load sites' data
    (
        zbar_2017,
        gamma_vals,
        z_2017,
        forestArea_2017_ha,
        theta_vals,
        site_theta_2017_df,
        site_gamma_2017_df,
    ) = load_site_data(site_num, norm_fac=norm_fac)

    num_sites = gamma_vals.size

    # Retrieving Stan data
    _, X_theta, N_theta, K_theta, G_theta, _ = theta_reg_data(
        num_sites, site_theta_2017_df
    )
    _, X_gamma, N_gamma, K_gamma, G_gamma = gamma_reg_data(
        num_sites, site_gamma_2017_df
    )

    # Save starting params
    uncertain_vals = np.concatenate((theta_vals, gamma_vals)).copy()
    uncertain_vals_old = np.concatenate((theta_vals, gamma_vals)).copy()

    # Retrieve z data for selected site(s)
    site_z_vals = z_2017

    # Collected Ensembles over all iterations; dictionary indexed by iteration number
    collected_ensembles = {}
    coe_ensembles = {}

    # Track error over iterations
    uncertain_vals_tracker = [uncertain_vals_old.copy()]
    abs_error_tracker = []
    percentage_error_tracker = []
    log_diff_error_tracker = []
    sol_val_X_tracker = []
    sol_val_Ua_tracker = []
    sol_val_Up_tracker = []
    sol_val_Um_tracker = []
    sol_val_Z_tracker = []
    sampling_time_tracker = []

    # Create dynamics matrices
    arr = np.cumsum(
        np.triu(np.ones((T, T))),
        axis=1,
    ).T
    Bdym = (1 - alpha) ** (arr - 1)
    Bdym[Bdym > 1] = 0.0
    Adym = np.arange(1, T + 1)
    alpha_p_Adym = np.power(1 - alpha, Adym)

    # Time step
    dt = T / N

    # Other placeholders!
    ds_vect = np.exp(-delta_t * np.arange(N + 1) * dt)
    ds_vect = np.reshape(ds_vect, (ds_vect.size, 1))

    # Results dictionary
    results = dict(
        num_sites=num_sites,
        tol=tol,
        T=T,
        N=N,
        norm_fac=norm_fac,
        delta_t=delta_t,
        alpha=alpha,
        kappa=kappa,
        pf=pf,
        pa=pa,
        xi=xi,
        zeta=zeta,
        sample_size=sample_size,
        final_sample_size=final_sample_size,
        weight=weight,
        output_dir=output_dir,
    )

    # Initialize error & iteration counter
    abs_error = np.infty
    percentage_error = np.infty
    cntr = 0

    # Loop until convergence
    while cntr < max_iter and percentage_error > tol:
        print(f"Optimization Iteration[{cntr+1}/{max_iter}]\n")

        # Flatten uncertain values
        uncertain_vals = np.asarray(uncertain_vals).flatten()

        # Unpacking uncertain values
        theta_vals = uncertain_vals[:num_sites].copy()
        gamma_vals = uncertain_vals[num_sites:].copy()

        print(f"Theta: {theta_vals}\n")
        print(f"Gamma: {gamma_vals}\n")

        # Computing carbon absorbed in start period
        x0_vals = gamma_vals * forestArea_2017_ha / norm_fac

        # Solve outer optimization problem
        (
            sol_val_X,
            sol_val_Up,
            sol_val_Um,
            sol_val_Z,
            sol_val_Ua,
        ) = solve_outer_optimization_problem(
            N=N,
            dt=dt,
            ds_vect=ds_vect,
            theta_vals=theta_vals,
            gamma_vals=gamma_vals,
            x0_vals=x0_vals,
            zbar_2017=zbar_2017,
            site_z_vals=site_z_vals,
            alpha=alpha,
            kappa=kappa,
            pf=pf,
            pa=pa,
            zeta=zeta,
        )

        # Update trackers
        sol_val_X_tracker.append(sol_val_X)
        sol_val_Ua_tracker.append(sol_val_Ua)
        sol_val_Up_tracker.append(sol_val_Up)
        sol_val_Um_tracker.append(sol_val_Um)
        sol_val_Z_tracker.append(sol_val_Z)

        # HMC sampling
        print("Starting HMC sampling...\n")
        model_data = dict(
            T=T,
            S=num_sites,
            norm_fac=norm_fac,
            alpha=alpha,
            sol_val_X=sol_val_X,
            sol_val_Ua=sol_val_Ua,
            sol_val_Up=sol_val_Up,
            zbar_2017=zbar_2017,
            forestArea_2017_ha=forestArea_2017_ha,
            alpha_p_Adym=alpha_p_Adym,
            Bdym=Bdym,
            ds_vect=ds_vect.flatten(),
            zeta=zeta,
            xi=xi,
            kappa=kappa,
            pa=pa,
            pf=pf,
            K_theta=K_theta,
            K_gamma=K_gamma,
            N_theta=N_theta,
            N_gamma=N_gamma,
            X_theta=X_theta,
            G_theta=G_theta,
            X_gamma=X_gamma,
            G_gamma=G_gamma,
            pa_2017=pa_2017,
            **_prior_hyperparams(num_sites, site_theta_2017_df, "theta"),
            **_prior_hyperparams(num_sites, site_gamma_2017_df, "gamma"),
        )

        # Compiling model
        sampler = stan.build(program_code=model_code, data=model_data, random_seed=1)
        print("Model compiled!\n")

        # Posterior sampling
        sampling_time = time.time()
        fit = sampler.sample(
            num_chains=num_chains, num_samples=sample_size, num_warmup=num_warmup
        )
        sampling_time = time.time() - sampling_time
        sampling_time_tracker.append(sampling_time)
        print(f"Finished sampling! Elapsed Time: {sampling_time} seconds\n")

        # Extract samples
        theta_post_samples = fit["theta"].T
        gamma_post_samples = fit["gamma"].T
        theta_coe_post_samples = fit["beta_theta"].T
        gamma_coe_post_samples = fit["beta_gamma"].T

        uncertainty_post_samples = np.concatenate(
            (theta_post_samples, gamma_post_samples), axis=1
        )

        uncertainty_coe_post_samples = np.concatenate(
            (theta_coe_post_samples, gamma_coe_post_samples), axis=1
        )

        # Update ensemble/tracker
        collected_ensembles.update({cntr: uncertainty_post_samples.copy()})
        coe_ensembles.update({cntr: uncertainty_coe_post_samples.copy()})

        print(f"Parameters from last iteration: {uncertain_vals_old}\n")
        print(
            f"""Parameters from current iteration:
            {np.mean(uncertainty_post_samples, axis=0)}\n"""
        )

        # Compute exponentially-smoothened new params
        uncertain_vals = (
            weight * np.mean(uncertainty_post_samples, axis=0)
            + (1 - weight) * uncertain_vals_old
        )

        uncertain_vals_tracker.append(uncertain_vals.copy())
        print(f"Updated uncertain values: {uncertain_vals}\n")

        # Evaluate error for convergence check
        # The percentage difference are changed to absolute difference
        abs_error = np.max(np.abs(uncertain_vals_old - uncertain_vals))
        percentage_error = np.max(
            np.abs(uncertain_vals_old - uncertain_vals) / uncertain_vals_old
        )

        abs_error_tracker.append(abs_error)
        percentage_error_tracker.append(percentage_error)

        print(
            f"""
            Iteration [{cntr+1:4d}]: Absolute Error = {abs_error},
            Percentage Error = {percentage_error}
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
                "percentage_error_tracker": np.asarray(percentage_error_tracker),
                "log_diff_error_tracker": np.asarray(log_diff_error_tracker),
                "uncertain_vals_tracker": np.asarray(uncertain_vals_tracker),
                "sampling_time_tracker": sampling_time_tracker,
                "collected_ensembles": collected_ensembles,
                "sol_val_X_tracker": sol_val_X_tracker,
                "sol_val_Ua_tracker": sol_val_Ua_tracker,
                "sol_val_Up_tracker": sol_val_Up_tracker,
                "sol_val_Um_tracker": sol_val_Um_tracker,
                "sol_val_Z_tracker": sol_val_Z_tracker,
                "coe_ensembles": coe_ensembles,
            }
        )

        # Save results (overwrite existing file)
        saveto = os.path.join(output_dir, "results.pcl")
        pickle.dump(results, open(saveto, "wb"))

    # Sample (densly) the final distribution
    print("Terminated. Sampling the final distribution...\n")
    fit = sampler.sample(num_chains=num_chains, num_samples=final_sample_size)

    theta_post_samples = fit["theta"].T
    gamma_post_samples = fit["gamma"].T
    theta_coe_post_samples = fit["beta_theta"].T
    gamma_coe_post_samples = fit["beta_gamma"].T

    final_sample = np.concatenate((theta_post_samples, gamma_post_samples), axis=1)
    final_sample_coe = np.concatenate(
        (theta_coe_post_samples, gamma_coe_post_samples), axis=1
    )

    results.update({"final_sample": final_sample})
    results.update({"final_sample_coe": final_sample_coe})

    # Save results (overwrite existing file)
    saveto = os.path.join(output_dir, "results.pcl")
    pickle.dump(results, open(saveto, "wb"))
    print(f"Results saved to {saveto}")

    return results


def _prior_hyperparams(num_sites, df, var):
    df = df.dropna()
    if var == "theta":
        # Get theta data
        y, X, _, _, _, W = theta_reg_data(num_sites, df)

        # Applying WLS weights
        y = W @ y
        X = W @ X

    elif var == "gamma":
        # Get gamma data
        y, X, _, _, _ = gamma_reg_data(num_sites, df)
    else:
        raise Exception("Argument `var` should be one of `theta`, `gamma`")

    inv_Lambda = np.linalg.inv(X.T @ X)
    mu = inv_Lambda @ X.T @ y
    a = (X.shape[0]) / 2
    b = 0.5 * (y.T @ y - mu.T @ X.T @ X @ mu)
    return {
        f"inv_Lambda_{var}": inv_Lambda,
        f"mu_{var}": mu,
        f"a_{var}": a,
        f"b_{var}": b,
    }
