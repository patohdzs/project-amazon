import time

import numpy as np
from cmdstanpy import CmdStanModel

from pysrc.sampling import baseline

from ..optimization import gams, gurobi
from ..sampling import gamma_adj_reg_data, theta_adj_reg_data
from ..sampling.baseline import baseline_hyperparams
from ..services.data_service import load_site_data
from ..services.file_service import stan_model_path


def sample(
    model_name,
    xi,
    pf,
    pa,
    weight,
    num_sites,
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
        gamma_vals,
        z_2017,
        forest_area_2017,
        theta_vals,
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
        pf=pf,
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

        solution = solve_planner_problem(
            T=T,
            theta=theta_vals,
            gamma=gamma_vals,
            x0=x0_vals,
            z0=z_2017,
            zbar=zbar_2017,
            dt=dt,
            pe=pf,
            pa=pa,
            alpha=alpha,
            delta=delta,
            kappa=kappa,
            zeta=zeta,
        )

        # Update trackers
        solution_tracker.append(solution)

        # HMC sampling
        print("Starting HMC sampling...\n")
        model_data = dict(
            T=T,
            S=num_sites,
            alpha=alpha,
            zbar_2017=zbar_2017,
            forest_area_2017=forest_area_2017,
            zeta=zeta,
            xi=xi,
            kappa=kappa,
            pa=pa,
            pa_2017=pa_2017,
            pf=pf,
            **solution,
            **_dynamics_matrices(T, dt, alpha, delta),
            **theta_adj_reg_data(num_sites, site_theta_df),
            **gamma_adj_reg_data(num_sites, site_gamma_df),
            **baseline_hyperparams(municipal_theta_df, "theta"),
            **baseline_hyperparams(municipal_gamma_df, "gamma"),
        )

        # Sampling from adjusted distribution
        sampling_time = time.time()
        fit = sampler.sample(
            data=model_data,
            **stan_kwargs,
        )
        sampling_time = time.time() - sampling_time
        print(f"Finished sampling! Elapsed Time: {sampling_time} seconds\n")
        print(fit.diagnose())

        # Update fit and sampling time trackers
        fit_tracker.append(fit.summary())
        sampling_time_tracker.append(sampling_time)

        # Extract samples
        theta_adj_samples = fit.stan_variable("theta")
        gamma_adj_samples = fit.stan_variable("gamma")
        theta_coe_adj_samples = fit.stan_variable("beta_theta")
        gamma_coe_adj_samples = fit.stan_variable("beta_gamma")

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
    stan_kwargs["iter_sampling"] = final_sample_size
    fit = sampler.sample(
        data=model_data,
        **stan_kwargs,
    )

    # Extract samples
    theta_adj_samples = fit.stan_variable("theta")
    gamma_adj_samples = fit.stan_variable("gamma")
    theta_coe_adj_samples = fit.stan_variable("beta_theta")
    gamma_coe_adj_samples = fit.stan_variable("beta_gamma")

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
