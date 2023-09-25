import os
import pickle

import numpy as np
import stan
from services.data_service import load_site_data
from services.file_service import stan_model_path
from solvers.casadi import solve_outer_optimization_problem


def sample_with_stan(
    model_name="basic_model.stan",
    # Configurations/Settings
    T=200,
    N=200,
    site_num=10,  # Number of sites(10, 25, 100, 1000)
    norm_fac=1e11,  # normalization factor consistent used in paper
    delta_t=0.02,
    alpha=0.045007414,
    kappa=2.094215255,
    pf=20.76,
    pa=44.75,
    xi=0.01,
    zeta=1.66e-4 * 1e11,  # use the same normalization factor
    max_iter=30,
    tol=0.001,
    sample_size=1000,
    final_sample_size=5_000,  # number of samples to collect after convergence
    num_chains=2,
    weight=0.25,  # <-- Not sure how this linear combination weighting helps!
    output_dir="output/stan_results",
):
    # Create the output directory
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    with open(stan_model_path(model_name)) as f:
        model_code = f.read()

    # Load sites' data
    (
        zbar_2017,
        gamma,
        z_2017,
        forestArea_2017_ha,
        theta,
        gamma_coe,
        gamma_coe_sd,
        theta_coe,
        theta_coe_sd,
        gamma_vcov_array,
        theta_vcov_array,
        site_theta_2017_df,
        site_gamma_2017_df,
    ) = load_site_data(
        site_num,
        norm_fac=norm_fac,
    )

    # Evaluate Gamma values
    gamma_vals = gamma
    num_sites = gamma.size

    # Theta Values
    theta_vals = theta
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
        print("Gamma: ", gamma_vals)
        print("Theta: ", theta_vals)

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

        sol_val_X_tracker.append(sol_val_X)
        sol_val_Ua_tracker.append(sol_val_Ua)
        sol_val_Up_tracker.append(sol_val_Up)
        sol_val_Um_tracker.append(sol_val_Um)
        sol_val_Z_tracker.append(sol_val_Z)

        print("Starting HMC sampling...")
        # Stan sampler
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
        )

        sampler = stan.build(program_code=model_code, data=model_data, random_seed=1)
        print("Model compiled!\n")

        fit = sampler.fixed_param(num_chains=num_chains, num_samples=sample_size)
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
            Iteration [{cntr+1:4d}]: Absolte Error = {abs_error},
            Percentage Error = {percentage_error}
            """
        )

        # Exchange gamma values (for future weighting/update & error evaluation)
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
    fit = sampler.fixed_param(num_chains=2, num_samples=final_sample_size)
    samples = fit.to_frame()

    # Final sample!!!
    final_sample = np.asarray(
        samples[[s for s in samples.columns if not s.endswith("_")]]
    )
    results.update({"final_sample": final_sample})

    # Save results (overwrite existing file)
    saveto = os.path.join(output_dir, "results.pcl")
    pickle.dump(results, open(saveto, "wb"))
    print(f"Results saved to {saveto}")

    return results
