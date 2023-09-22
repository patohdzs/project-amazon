import os
import pickle
import time
from functools import partial

import casadi
import matplotlib.pyplot as plt
import numpy as np
from mcmc.hmc import create_hmc_sampler
from services.data_service import load_site_data
from solvers import _DEBUG, construct_block_matrix, log_density_function
from utils.text import decorate_text


def solve_with_casadi(
    # Configurations/Settings
    site_num=25,  # Number of sites(10, 25, 100, 1000)
    norm_fac=1e11,  # normalization factor consistent used in paper
    delta_t=0.02,
    alpha=0.045007414,
    kappa=2.094215255,
    pf=20.76,
    pa=44.75,
    xi=0.01,
    zeta=1.66e-4 * 1e11,  # use the same normalization factor
    max_iter=20000,
    tol=0.001,
    T=200,
    N=200,
    sample_size=1000,
    mode_as_solution=False,  # If true, use the modeas solution for gamma
    final_sample_size=5_000,  # number of samples to collect after convergence
    two_param_uncertainty=True,
    weight=0.25,  # <-- Not sure how this linear combination weighting helps!
    output_dir="Casadi_Results",
    mix_in=2,
    symplectic_integrator_num_steps=10,
    mass_matrix_weight=0.1,
    stepsize=0.1,
    scale=0.0,
):
    """
    Main function to solve the bilievel optimization problem using casadi for
    the outer optimization problem.

    :param float tol: convergence tolerance
    :param T:
    :param N:
    """

    # Create the output directory
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

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

    # Print the data
    print("data loaded")

    site_theta_2017_df = site_theta_2017_df
    site_gamma_2017_df = site_gamma_2017_df

    # Evaluate Gamma values
    gamma_vals = gamma
    size = gamma.size

    # Theta Values
    theta_vals = theta

    # Retrieve z data for selected site(s)
    site_z_vals = z_2017

    # Coes value
    gamma_coe_vals = gamma_coe
    theta_coe_vals = theta_coe

    # Construct the block matrix
    block_matrix = construct_block_matrix(theta_vcov_array, gamma_vcov_array)

    # Two parameter uncertainty (both theta and gamma)
    vals = np.concatenate((theta_coe_vals, gamma_coe_vals))
    size_coe = vals.size

    # Initialize Gamma Values
    uncertain_vals = vals.copy()
    uncertain_vals_mean = vals.copy()
    uncertain_vals_old = vals.copy()

    # Householder to track sampled gamma values
    uncertain_vals_tracker = [uncertain_vals.copy()]
    uncertain_SD_tracker = [block_matrix.copy()]

    mass_matrix = 10000 * np.diag(1 / np.concatenate((theta_coe_sd, gamma_coe_sd)) ** 2)
    mass_matrix_tracker = [mass_matrix.copy()]

    # Collected Ensembles over all iterations; dictionary indexed by iteration number
    collected_ensembles = {}

    # Track error over iterations
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

    # Initialize Blocks of the A matrix those won't change
    A = np.zeros((size + 2, size + 2))
    Ax = np.zeros(size + 2)

    # Construct Matrix B
    B = np.eye(N=size + 2, M=size, k=0)
    B = casadi.sparsify(B)

    # Construct Matrxi D constant blocks
    D = np.zeros((size + 2, size))

    # time step!
    dt = T / N

    # Other placeholders!
    ds_vect = np.exp(-delta_t * np.arange(N + 1) * dt)
    ds_vect = np.reshape(ds_vect, (ds_vect.size, 1))

    # Results dictionary
    results = dict(
        size=size,
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
        mode_as_solution=mode_as_solution,
        mix_in=mix_in,
        symplectic_integrator_num_steps=symplectic_integrator_num_steps,
        two_param_uncertainty=two_param_uncertainty,
        weight=weight,
        mass_matrix_weight=mass_matrix_weight,
        output_dir=output_dir,
        stepsize=stepsize,
    )

    # Initialize error & iteration counter
    abs_error = np.infty
    percentage_error = np.infty
    cntr = 0

    # Loop until convergence
    while cntr < max_iter and percentage_error > tol:
        print(decorate_text(f"Optimization Iteration[{cntr+1}/{max_iter}]"))

        print("uncertain_vals used in optimation: ", uncertain_vals)

        x0_vals = gamma_vals * forestArea_2017_ha / norm_fac

        # Construct Matrix A from new uncertain_vals
        A[:-2, :] = 0.0
        Ax[0:size] = -alpha * gamma_vals
        Ax[-1] = alpha * np.sum(gamma_vals * zbar_2017)
        Ax[-2] = -alpha
        A[-2, :] = Ax
        A[-1, :] = 0.0
        A = casadi.sparsify(A)

        # Construct Matrix D from new uncertain_vals
        D[:, :] = 0.0
        D[-2, :] = -gamma_vals
        D = casadi.sparsify(D)

        # Define the right hand side (symbolic here) as a function of gamma
        gamma = casadi.MX.sym("gamma", size + 2)
        up = casadi.MX.sym("up", size)
        um = casadi.MX.sym("um", size)

        rhs = (A @ gamma + B @ (up - um) + D @ up) * dt + gamma
        f = casadi.Function("f", [gamma, um, up], [rhs])

        ## Define an optimizer and initialize it, and set constraints
        opti = casadi.Opti()

        # Decision variables for states
        X = opti.variable(size + 2, N + 1)

        # Aliases for states
        Up = opti.variable(size, N)
        Um = opti.variable(size, N)
        Ua = opti.variable(1, N)

        # 1.2: Parameter for initial state
        ic = opti.parameter(size + 2)

        # Gap-closing shooting constraints
        for k in range(N):
            opti.subject_to(X[:, k + 1] == f(X[:, k], Um[:, k], Up[:, k]))

        # Initial and terminal constraints
        opti.subject_to(X[:, 0] == ic)
        opti.subject_to(opti.bounded(0, X[0:size, :], zbar_2017[0:size]))

        # Objective: regularization of controls
        for k in range(size):
            opti.subject_to(opti.bounded(0, Um[k, :], casadi.inf))
            opti.subject_to(opti.bounded(0, Up[k, :], casadi.inf))

        opti.subject_to(Ua == casadi.sum1(Up + Um) ** 2)

        term1 = casadi.sum2(ds_vect[0:N, :].T * Ua * zeta / 2)
        term2 = -casadi.sum2(ds_vect[0:N, :].T * (pf * (X[-2, 1:] - X[-2, 0:-1])))
        term3 = -casadi.sum2(
            ds_vect.T * casadi.sum1((pa * theta_vals - pf * kappa) * X[0:size, :])
        )

        opti.minimize(term1 + term2 + term3)

        # Solve optimization problem
        options = dict()
        options["print_time"] = True
        options["expand"] = True
        options["ipopt"] = {
            "print_level": 1,
            "fast_step_computation": "yes",
            "mu_allow_fast_monotone_decrease": "yes",
            "warm_start_init_point": "yes",
        }
        opti.solver("ipopt", options)

        opti.set_value(
            ic,
            casadi.vertcat(site_z_vals, np.sum(x0_vals), 1),
        )

        if _DEBUG:
            print("ic: ", ic)
            print("site_z_vals: ", site_z_vals)
            print("x0_vals: ", x0_vals)
            print(
                "casadi.vertcat(site_z_vals,np.sum(x0_vals),1): ",
                casadi.vertcat(site_z_vals, np.sum(x0_vals), 1),
            )

        # TODO: Discuss with Daniel how this is taking too long, not the sampling!
        print("Solving the outer optimization problem...")
        start_time = time.time()
        sol = opti.solve()
        print(f"Done; time taken {time.time()-start_time} seconds...")

        print("sol.value(X)", sol.value(X))
        print("sol.value(Ua)", sol.value(Ua))
        print("sol.value(Up)", sol.value(Up))
        print("sol.value(Um)", sol.value(Um))

        # Extract information from the solver
        N = X.shape[1] - 1
        sol_val_X = sol.value(X)
        sol_val_Up = sol.value(Up)
        sol_val_Um = sol.value(Um)
        sol_val_Z = sol_val_Up - sol_val_Um
        sol_val_Ua = sol.value(Ua)

        sol_val_X_tracker.append(sol_val_X)
        sol_val_Ua_tracker.append(sol_val_Ua)
        sol_val_Up_tracker.append(sol_val_Up)
        sol_val_Um_tracker.append(sol_val_Um)
        sol_val_Z_tracker.append(sol_val_Z)

        log_density = partial(
            log_density_function,
            uncertain_vals_mean=uncertain_vals_mean,
            block_matrix=block_matrix,
            alpha=alpha,
            N=N,
            sol_val_X=sol_val_X,
            sol_val_Ua=sol_val_Ua,
            sol_val_Up=sol_val_Up,
            zbar_2017=zbar_2017,
            forestArea_2017_ha=forestArea_2017_ha,
            norm_fac=norm_fac,
            alpha_p_Adym=alpha_p_Adym,
            Bdym=Bdym,
            leng=leng,
            T=T,
            ds_vect=ds_vect,
            zeta=zeta,
            xi=xi,
            kappa=kappa,
            pa=pa,
            pf=pf,
            site_theta_2017_df=site_theta_2017_df,
            site_gamma_2017_df=site_gamma_2017_df,
        )

        print("Starting HMC sampling...")
        sampler = create_hmc_sampler(
            size=size_coe,
            log_density=log_density,
            burn_in=100,
            mix_in=mix_in,
            symplectic_integrator="verlet",
            symplectic_integrator_stepsize=stepsize,
            symplectic_integrator_num_steps=symplectic_integrator_num_steps,
            constraint_test=lambda x: True if np.max(x >= 0) else False,
            mass_matrix=mass_matrix,
        )

        sampling_results = sampler.start_MCMC_sampling(
            sample_size=sample_size,
            initial_state=uncertain_vals,
            verbose=True,
        )

        uncertainty_post_samples = np.asarray(sampling_results["collected_ensemble"])
        uncertainty_map_estimate = sampling_results["map_estimate"]

        # Update ensemble/tracker
        collected_ensembles.update({cntr: uncertainty_post_samples.copy()})

        if mode_as_solution:
            print("uncertain values from last iteration: ", uncertain_vals_old)
            print("uncertain values from this iteration: ", uncertainty_map_estimate)
            if scale == 1.0:
                uncertainty_map_estimate[:size] *= thetaSD  # noqa: F405
                uncertainty_map_estimate[size:] *= gammaSD  # noqa: F405
                print(
                    "scale back uncertain values from this iteration: ",
                    uncertainty_map_estimate,
                )

            uncertain_vals = (
                weight * uncertainty_map_estimate + (1 - weight) * uncertain_vals_old
            )
            print("updated uncertain values: ", uncertain_vals)

        else:
            print("uncertain values from last iteration: ", uncertain_vals_old)
            print(
                "uncertain values from this iteration: ",
                np.mean(uncertainty_post_samples, axis=0),
            )

            uncertain_vals = (
                weight * np.mean(uncertainty_post_samples, axis=0)
                + (1 - weight) * uncertain_vals_old
            )

            print("updated uncertain values: ", uncertain_vals)

        uncertain_vals_tracker.append(uncertain_vals.copy())

        ## to do
        theta_coe_subset = uncertainty_post_samples[:, :8]
        gamma_coe_subset = uncertainty_post_samples[:, 8:]
        theta_vcov_array = np.cov(theta_coe_subset, rowvar=False)
        gamma_vcov_array = np.cov(gamma_coe_subset, rowvar=False)

        # Construct the block matrix
        block_matrix_post = construct_block_matrix(theta_vcov_array, gamma_vcov_array)

        uncertain_post_SD = block_matrix_post
        print(
            "uncertain value standard deviation from this iteration: ",
            uncertain_post_SD,
        )
        uncertain_SD_tracker.append(uncertain_post_SD.copy())
        print("updated mass matrix:", mass_matrix)

        mass_matrix_tracker.append(mass_matrix.copy())

        # Evaluate error for convergence check
        # The percentage difference are changed to absolute difference
        abs_error = np.max(np.abs(uncertain_vals_old - uncertain_vals))
        percentage_error = np.max(
            np.abs(uncertain_vals_old - uncertain_vals) / uncertain_vals_old
        )

        abs_error_tracker.append(abs_error)
        percentage_error_tracker.append(percentage_error)

        print(
            decorate_text(
                f"""
                Iteration [{cntr+1:4d}]: Absolte Error = {abs_error},
                Percentage Error = {percentage_error}
                """
            )
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

        # Extensive plotting for monitoring; not needed really!
        if False:
            plt.plot(uncertain_vals_tracker[-2], label=r"Old $\gamma$")
            plt.plot(uncertain_vals_tracker[-1], label=r"New $\gamma$")
            plt.legend()
            plt.show()

            for j in range(size):
                plt.hist(uncertainty_post_samples[:, j], bins=50)
                plt.title(f"Iteration {cntr}; Site {j+1}")
                plt.show()

    print("Terminated. Sampling the final distribution...")
    # Sample (densly) the final distribution
    final_sample = sampler.sample(
        sample_size=final_sample_size,
        initial_state=uncertain_vals,
        verbose=True,
    )
    final_sample = np.asarray(final_sample)
    results.update({"final_sample": final_sample})

    # Save results (overwrite existing file)
    saveto = os.path.join(output_dir, "results.pcl")
    pickle.dump(results, open(saveto, "wb"))
    print(f"Results saved to {saveto}")

    return results
