import os
import pickle

import numpy as np
import stan
from services.data_service import load_site_data
from services.file_service import stan_model_path
from solvers.casadi import solve_outer_optimization_problem


def theta_reg_data(num_sites, theta_df):
    theta_df = theta_df[theta_df["zbar_2017_muni"].notna()]

    X_theta = theta_df.iloc[:, 1:9].to_numpy()
    N_theta, K_theta = X_theta.shape

    G_theta = np.array(
        [(theta_df["id"].to_numpy() == i).astype(int) for i in range(1, num_sites + 1)]
    )
    G_theta = theta_df["zbar_2017_muni"].to_numpy() * G_theta
    G_theta = G_theta / G_theta.sum(axis=1, keepdims=True)

    return X_theta, N_theta, K_theta, G_theta


def gamma_reg_data(num_sites, gamma_df):
    X_gamma = gamma_df.iloc[:, 1:6].to_numpy()
    N_gamma, K_gamma = X_gamma.shape
    G_gamma = np.array(
        [(gamma_df["id"].to_numpy() == i).astype(int) for i in range(1, num_sites + 1)]
    )
    G_gamma = G_gamma / G_gamma.sum(axis=1, keepdims=True)

    return X_gamma, N_gamma, K_gamma, G_gamma


def sample_with_stan(
    model_name,
    output_dir,
    xi,
    pf,
    pa,
    weight,
    site_num,
    T,
    N=200,
    norm_fac=1e11,
    delta_t=0.02,
    alpha=0.045007414,
    kappa=2.094215255,
    zeta=1.66e-4 * 1e11,  # use the same normalization factor
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
    with open(stan_model_path(model_name)) as f:
        model_code = f.read()

    # Load sites' data
    (
        zbar_2017,
        gamma_vals,
        z_2017,
        forestArea_2017_ha,
        theta_vals,
        gamma_coe,
        gamma_coe_sd,
        theta_coe,
        theta_coe_sd,
        gamma_vcov_array,
        theta_vcov_array,
        site_theta_2017_df,
        site_gamma_2017_df,
    ) = load_site_data(site_num, norm_fac=norm_fac)

    num_sites = gamma_vals.size

    # Splitting data
    X_theta, N_theta, K_theta, G_theta = theta_reg_data(num_sites, site_theta_2017_df)
    X_gamma, N_gamma, K_gamma, G_gamma = gamma_reg_data(num_sites, site_gamma_2017_df)

    # Save starting params
    uncertain_vals = np.concatenate((theta_vals, gamma_vals)).copy()
    uncertain_vals_old = np.concatenate((theta_vals, gamma_vals)).copy()

    # Retrieve z data for selected site(s)
    site_z_vals = z_2017

    # Collected Ensembles over all iterations; dictionary indexed by iteration number
    collected_ensembles = {}

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

    # Update this parameter (leng) once figured out where it is coming from
    leng = 200
    arr = np.cumsum(
        np.triu(np.ones((leng, leng))),
        axis=1,
    ).T
    Bdym = (1 - alpha) ** (arr - 1)
    Bdym[Bdym > 1] = 0.0
    Adym = np.arange(1, leng + 1)
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
        beta_theta_prior_mean=theta_coe,
        beta_theta_prior_vcov=theta_vcov_array,
        beta_gamma_prior_mean=gamma_coe,
        beta_gamma_prior_vcov=gamma_vcov_array,
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
        print(f"Optimization Iteration[{cntr+1}/{max_iter}]")

        # Flatten uncertain values
        uncertain_vals = np.asarray(uncertain_vals).flatten()

        # Unpacking uncertain values
        theta_vals = uncertain_vals[:num_sites].copy()
        gamma_vals = uncertain_vals[num_sites:].copy()

        print("Theta: ", theta_vals)
        print("Gamma: ", gamma_vals)

        # Unpacking uncertain values
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
        print("Starting HMC sampling...")
        model_data = dict(
            T=T,
            S=num_sites,
            K_theta=K_theta,
            K_gamma=K_gamma,
            N_theta=N_theta,
            N_gamma=N_gamma,
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
            X_theta=X_theta,
            G_theta=G_theta,
            X_gamma=X_gamma,
            G_gamma=G_gamma,
            beta_theta_prior_mean=theta_coe,
            beta_theta_prior_vcov=theta_vcov_array,
            beta_gamma_prior_mean=gamma_coe,
            beta_gamma_prior_vcov=gamma_vcov_array,
        )

        # Compiling model
        sampler = stan.build(program_code=model_code, data=model_data, random_seed=1)
        print("Model compiled!\n")

        # Posterior sampling
        fit = sampler.sample(
            num_chains=num_chains, num_samples=sample_size, num_warmup=num_warmup
        )
        print("Finished sampling!\n")

        samples = fit.to_frame()

        theta_post_samples = np.asarray(
            samples[[s for s in samples.columns if s.startswith("theta")]]
        )

        gamma_post_samples = np.asarray(
            samples[[s for s in samples.columns if s.startswith("gamma")]]
        )

        uncertainty_post_samples = np.concatenate(
            (theta_post_samples, gamma_post_samples), axis=1
        )

        # Update ensemble/tracker
        collected_ensembles.update({cntr: uncertainty_post_samples.copy()})

        print("Parameters from last iteration: ", uncertain_vals_old)
        print(
            "Parameters from this iteration: ",
            np.mean(uncertainty_post_samples, axis=0),
        )

        uncertain_vals = (
            weight * np.mean(uncertainty_post_samples, axis=0)
            + (1 - weight) * uncertain_vals_old
        )

        print("Updated uncertain values: ", uncertain_vals)

        uncertain_vals_tracker.append(uncertain_vals.copy())
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
                "collected_ensembles": collected_ensembles,
                "sol_val_X_tracker": sol_val_X_tracker,
                "sol_val_Ua_tracker": sol_val_Ua_tracker,
                "sol_val_Up_tracker": sol_val_Up_tracker,
                "sol_val_Um_tracker": sol_val_Um_tracker,
                "sol_val_Z_tracker": sol_val_Z_tracker,
            }
        )

        # Save results (overwrite existing file)
        saveto = os.path.join(output_dir, "results.pcl")
        pickle.dump(results, open(saveto, "wb"))

    # Sample (densly) the final distribution
    print("Terminated. Sampling the final distribution...")
    fit = sampler.sample(num_chains=num_chains, num_samples=final_sample_size)
    samples = fit.to_frame()

    # Final sample!!!
    theta_post_samples = np.asarray(
        samples[[s for s in samples.columns if s.startswith("theta")]]
    )
    gamma_post_samples = np.asarray(
        samples[[s for s in samples.columns if s.startswith("gamma")]]
    )
    final_sample = np.concatenate((theta_post_samples, gamma_post_samples), axis=1)

    results.update({"final_sample": final_sample})

    # Save results (overwrite existing file)
    saveto = os.path.join(output_dir, "results.pcl")
    pickle.dump(results, open(saveto, "wb"))
    print(f"Results saved to {saveto}")

    return results
