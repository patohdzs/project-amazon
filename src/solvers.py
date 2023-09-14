#!/usr/bin/env python
# coding: utf-8

"""
This module presents solution approaches for solving the bilevel optimization problem with
MCMC sampling for the inner optimization problem.
The outer optimization problem is either solved with CASADI or GAMS.
"""


# Import Required Packages
# ========================
import os, sys
import pickle
import time
import shutil
import getpass
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy.ma as ma
from scipy.stats import truncnorm

# MCMC (HMC) sampling routines

mcmc_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "src/mcmc")
if mcmc_path not in sys.path:
    sys.path.append(mcmc_path)

# Data Hanlder (.data_handlers.load_site_data)
sys.path.append(os.path.abspath("src"))
from services.data_service import load_site_data


# check if casadi is available; delay exception raise to the call
try:
    import casadi
except ImportError:
    casadi = None

# check if gams is available; delay exception raise to the call
try:
    from gams import GamsWorkspace

    if GamsWorkspace.api_major_rel_number < 42:  # old API structure
        import gdxcc as gdx
        from gams import *
        import gamstransfer as gt
    else:  # new API structure
        import gams.core.gdx as gdx
        from gams.control import *
        import gams.transfer as gt
except ImportError:
    gams = None


# Usefule Variables:
_DEBUG = False
_GAMS_SYSTEM_LOADER = os.path.join(
    os.path.dirname(__file__), "_gams_system_directory.dat"
)


class TextColor:
    """Class with text colors and font settings for printing"""

    PURPLE = "\033[95m"
    CYAN = "\033[96m"
    DARKCYAN = "\033[36m"
    BLUE = "\033[94m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    RED = "\033[91m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"
    END = "\033[0m"


def decorate_text(text):
    sep = "*" * 48
    return f"{TextColor.BOLD}{TextColor.DARKCYAN}\n\n{sep}\n\t{text}\n{sep}\n{TextColor.END}\n"


def get_gams_system_directory(
    filepath=_GAMS_SYSTEM_LOADER,
):
    """
    Load the GAMS system directory from file or ask the user about it

    An example is that the GAMS System Directory could be (this is what Was hard coded):
    `/Library/Frameworks/GAMS.framework/Versions/43/Resources`
    """
    if not os.path.isfile(filepath):
        gams_sys_dir = None

    else:
        # Load the path from file and validate
        with open(filepath, "r") as f_id:
            gams_sys_dir = f_id.read().strip(" \n")

        if not gams_sys_dir:
            # The filepath is empty and will be overwritten
            # os.remove(filepath)
            gams_sys_dir = None

        elif not os.path.isdir(gams_sys_dir):
            # The file contains a path to an invalid directory
            print(
                f"The GAMS system directory below is not valid\n"
                f"Invalid GAMS System Dir: '{gams_sys_dir}'"
            )
            gams_sys_dir = None

    if gams_sys_dir is None:
        # Either file does not exist or path in it is invalid.
        # Ask user for a valid path, then write it to file and validate
        prompt = f"\n**\nPlease input FULL path to GAMS system directory/resources.\n"
        prompt += f"For example: '/Library/Frameworks/GAMS.framework/Versions/43/Resources'\n$ "
        gams_sys_dir = input(prompt).strip(""" \n" '  """)
        # Write it to file and recurse
        with open(filepath, "w") as f_id:
            f_id.write(gams_sys_dir)

        # Recurse to validate the path
        return get_gams_system_directory(
            filepath=filepath,
        )

    return gams_sys_dir


def log_density_function(
    uncertain_vals,
    uncertain_vals_mean,
    theta_vals,
    N,
    site_precisions,
    alpha,
    sol_val_X,
    sol_val_Ua,
    sol_val_Up,
    zbar_2017,
    forestArea_2017_ha,
    norm_fac,
    alpha_p_Adym,
    Bdym,
    leng,
    T,
    ds_vect,
    zeta,
    xi,
    kappa,
    pa,
    pf,
    two_param_uncertainty=True,
    thetaSD=None,
    gammaSD=None,
    scale=0.0,
):
    """
    Define a function to evaluate log-density of the objective/posterior distribution
    Some of the input parameters are updated at each cycle of the outer loop (optimization loop),
    and it becomes then easier/cheaper to udpate the function stamp and keep it separate here

    Note that the log-density is the logarithm of the target density discarding any normalization factor
    """
    if not two_param_uncertainty:
        # One parameter (gamma) uncertainty
        ds_vect = np.asarray(ds_vect).flatten()
        uncertain_vals = np.asarray(uncertain_vals).flatten()
        size = theta_vals.size
        x0_vals = uncertain_vals.T.dot(forestArea_2017_ha) / norm_fac
        X_zero = np.sum(x0_vals) * np.ones(leng)

        # shifted_X = zbar_2017 - sol.value(X)[0:size, :-1]
        shifted_X = sol_val_X[0:size, :-1].copy()
        for j in range(N):
            shifted_X[:, j] = zbar_2017 - shifted_X[:, j]
        omega = np.dot(uncertain_vals, alpha * shifted_X - sol_val_Up)

        X_dym = np.zeros(T + 1)
        X_dym[0] = np.sum(x0_vals)
        X_dym[1:] = alpha_p_Adym * X_zero + np.dot(Bdym, omega.T)

        z_shifted_X = sol_val_X[0:size, :].copy()
        scl = pa * theta_vals - pf * kappa
        for j in range(N + 1):
            z_shifted_X[:, j] *= scl

        term_1 = -np.sum(ds_vect[0:T] * sol_val_Ua) * zeta / 2
        term_2 = np.sum(ds_vect[0:T] * (X_dym[1:] - X_dym[0:-1])) * pf
        term_3 = np.sum(ds_vect * np.sum(z_shifted_X, axis=0))

        obj_val = term_1 + term_2 + term_3

        uncertain_vals_dev = uncertain_vals - uncertain_vals_mean
        norm_log_prob = -0.5 * np.dot(
            uncertain_vals_dev, site_precisions.dot(uncertain_vals_dev)
        )
        log_density_val = -1.0 / xi * obj_val + norm_log_prob

        log_density_val = float(log_density_val)

        if _DEBUG:
            print("Term 1: ", term_1)
            print("Term 2: ", term_2)
            print("Term 3: ", term_3)
            print("obj_val: ", obj_val)
            print("norm_log_prob", norm_log_prob)
            print("log_density_val", log_density_val)
    else:
        # Two parameter uncertainty (both theta and gamma)
        ds_vect = np.asarray(ds_vect).flatten()
        uncertain_vals = np.asarray(uncertain_vals).flatten()
        size = theta_vals.size
        if scale == 1.0:
            uncertain_vals[:size] = uncertain_vals[:size] * thetaSD
            uncertain_vals[size:] = uncertain_vals[size:] * gammaSD
        # print("uncertain_vals used in log density:", uncertain_vals)
        x0_vals = uncertain_vals[size:].T.dot(forestArea_2017_ha) / norm_fac
        X_zero = np.sum(x0_vals) * np.ones(leng)

        # shifted_X = zbar_2017 - sol.value(X)[0:size, :-1]
        shifted_X = sol_val_X[0:size, :-1].copy()
        for j in range(N):
            shifted_X[:, j] = zbar_2017 - shifted_X[:, j]
        omega = np.dot(uncertain_vals[size:], alpha * shifted_X - sol_val_Up)

        X_dym = np.zeros(T + 1)
        X_dym[0] = np.sum(x0_vals)
        X_dym[1:] = alpha_p_Adym * X_zero + np.dot(Bdym, omega.T)

        z_shifted_X = sol_val_X[0:size, :].copy()
        scl = pa * uncertain_vals[:size] - pf * kappa
        for j in range(N + 1):
            z_shifted_X[:, j] *= scl

        term_1 = -np.sum(ds_vect[0:T] * sol_val_Ua) * zeta / 2
        term_2 = np.sum(ds_vect[0:T] * (X_dym[1:] - X_dym[0:-1])) * pf
        term_3 = np.sum(ds_vect * np.sum(z_shifted_X, axis=0))

        obj_val = term_1 + term_2 + term_3
        if scale == 1.0:
            uncertain_vals[:size] = uncertain_vals[:size] / thetaSD
            uncertain_vals[size:] = uncertain_vals[size:] / gammaSD
        uncertain_vals_dev = uncertain_vals - uncertain_vals_mean
        norm_log_prob = -0.5 * np.dot(
            uncertain_vals_dev, site_precisions.dot(uncertain_vals_dev)
        )
        log_density_val = -1.0 / xi * obj_val + norm_log_prob
        log_density_val = float(log_density_val)

    return log_density_val


def solve_with_casadi(
    # Configurations/Settings
    site_num=25,  # Number of sites(10, 25, 100, 1000)
    norm_fac=1e9,
    delta_t=0.02,
    alpha=0.045007414,
    kappa=2.094215255,
    pf=20.76,
    pa=44.75,
    xi=0.01,
    zeta=1.66e-4 * 1e9,  # zeta := 1.66e-4*norm_fac  #
    #
    max_iter=20000,
    tol=0.001,
    T=200,
    N=200,
    #
    sample_size=1000,  # simulations before convergence (to evaluate the mean)
    mode_as_solution=False,  # If true, use the mode (point of high probability) as solution for gamma
    final_sample_size=25_000,  # number of samples to collect after convergence
    two_param_uncertainty=True,
    weight=0.25,  # <-- Not sure how this linear combination weighting helps!
    output_dir="Casadi_Results",
    mix_in=2,
    mass_matrix_theta_scale=1.0,
    mass_matrix_gamma_scale=1.0,
    symplectic_integrator_num_steps=10,
    mass_matrix_weight=0.1,
    stepsize=0.1,
    scale=0.0,
    mode=1.0,
):
    if (mode == 1.0) or (mode == 2.0):
        from mcmc_sampling_original import create_hmc_sampler
    elif mode == 3.0:
        from mcmc_sampling_reflection import create_hmc_sampler
    """
    Main function to solve the bilievel optimization problem using casadi for the outer optimization problem.

    :param float tol: convergence tolerance
    :param T:
    :param N:
    """

    # if casadi is None:
    #     print("Failed to import CASADI. This function requires CASADI to be installed! ")
    #     raise ImportError

    # Create the output directory
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # Load sites' data
    (
        zbar_2017,
        gamma,
        gammaSD,
        z_2017,
        forestArea_2017_ha,
        theta,
        thetaSD,
    ) = load_site_data(
        site_num,
        norm_fac=norm_fac,
    )

    # Print the data
    print("data loaded")
    print("theta", theta, "thetaSD", thetaSD)
    print("gamma", gamma, "gammaSD", gammaSD)

    # Evaluate Gamma values ()
    gamma_1_vals = gamma - gammaSD
    gamma_2_vals = gamma + gammaSD
    gamma_vals = gamma
    size = gamma.size

    # Theta Values
    theta_vals = theta
    tot_size = gamma.size + theta.size

    # Retrieve z data for selected site(s)
    site_z_vals = z_2017

    if not two_param_uncertainty:
        # One parameter (gamma) uncertainty
        # Evaluate mean and covariances from site data
        site_stdev = gammaSD
        site_covariances = np.diag(np.power(site_stdev, 2))
        site_precisions = np.linalg.inv(site_covariances)
        site_mean = gamma_1_vals / 2 + gamma_2_vals / 2

        # Initialize Gamma Values
        uncertain_vals = gamma.copy()
        uncertain_vals_mean = gamma.copy()
        uncertain_vals_old = gamma.copy()

    else:
        # Two parameter uncertainty (both theta and gamma)
        vals = np.concatenate((theta_vals, gamma_vals))
        # Evaluate mean and covariances from site data
        site_stdev = np.concatenate((thetaSD, gammaSD))
        site_covariances = np.diag(np.power(site_stdev, 2))
        site_precisions = np.linalg.inv(site_covariances)
        site_mean = vals

        # Initialize Gamma Values
        uncertain_vals = vals.copy()
        uncertain_vals_mean = vals.copy()
        uncertain_vals_old = vals.copy()

    # Householder to track sampled gamma values
    # uncertain_vals_tracker       = np.empty((uncertain_vals.size, sample_size+1))
    # uncertain_vals_tracker[:, 0] = uncertain_vals.copy()
    uncertain_vals_tracker = [uncertain_vals.copy()]
    uncertain_SD_tracker = [np.concatenate((thetaSD, gammaSD)).copy()]

    if scale == 0.0:
        mass_matrix_theta = 1 / (thetaSD**2)
        mass_matrix_gamma = 1 / (gammaSD**2)
    else:
        mass_matrix_theta = np.ones(thetaSD.shape)
        mass_matrix_gamma = np.ones(gammaSD.shape)
    mass_matrix = np.concatenate((mass_matrix_theta, mass_matrix_gamma))
    mass_matrix_tracker = [mass_matrix.copy()]
    print("mass_matrix initialization:", mass_matrix)

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
    count_Ratio_tracker = []
    site_count_Ratio_tracker = []
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
        mass_matrix_theta_scale=mass_matrix_theta_scale,
        mass_matrix_gamma_scale=mass_matrix_gamma_scale,
        symplectic_integrator_num_steps=symplectic_integrator_num_steps,
        two_param_uncertainty=two_param_uncertainty,
        gammaSD=gammaSD,
        thetaSD=thetaSD,
        weight=weight,
        mass_matrix_weight=mass_matrix_weight,
        output_dir=output_dir,
        stepsize=stepsize,
    )

    # Initialize error & iteration counter
    abs_error = np.infty
    percentage_error = np.infty
    log_diff_error = np.infty
    cntr = 0

    # Loop until convergence
    while cntr < max_iter and percentage_error > tol:
        print(decorate_text(f"Optimization Iteration[{cntr+1}/{max_iter}]"))

        print("uncertain_vals used in optimation: ", uncertain_vals)
        if not two_param_uncertainty:
            # One parameter (gamma) uncertainty
            # Update x0
            x0_vals = uncertain_vals * forestArea_2017_ha / norm_fac
            # Construct Matrix A from new uncertain_vals
            A[:-2, :] = 0.0
            Ax[0:size] = -alpha * uncertain_vals[0:size]
            Ax[-1] = alpha * np.sum(uncertain_vals * zbar_2017)
            Ax[-2] = -alpha
            A[-2, :] = Ax
            A[-1, :] = 0.0
            A = casadi.sparsify(A)

            # Construct Matrix D from new uncertain_vals
            D[:, :] = 0.0
            D[-2, :] = -uncertain_vals
            D = casadi.sparsify(D)

        else:
            # Two parameter uncertainty (both theta and gamma)

            x0_vals = uncertain_vals[size:] * forestArea_2017_ha / norm_fac

            # Construct Matrix A from new uncertain_vals
            A[:-2, :] = 0.0
            Ax[0:size] = -alpha * uncertain_vals[size:]
            Ax[-1] = alpha * np.sum(uncertain_vals[size:] * zbar_2017)
            Ax[-2] = -alpha
            A[-2, :] = Ax
            A[-1, :] = 0.0
            A = casadi.sparsify(A)

            # Construct Matrix D from new uncertain_vals
            D[:, :] = 0.0
            D[-2, :] = -uncertain_vals[size:]
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

        if not two_param_uncertainty:
            # One parameter (gamma) uncertainty

            # Set teh optimization problem
            term1 = casadi.sum2(ds_vect[0:N, :].T * Ua * zeta / 2)
            term2 = -casadi.sum2(ds_vect[0:N, :].T * (pf * (X[-2, 1:] - X[-2, 0:-1])))
            term3 = -casadi.sum2(
                ds_vect.T * casadi.sum1((pa * theta_vals - pf * kappa) * X[0:size, :])
            )

        else:
            # Two parameter uncertainty (both theta and gamma)

            term1 = casadi.sum2(ds_vect[0:N, :].T * Ua * zeta / 2)
            term2 = -casadi.sum2(ds_vect[0:N, :].T * (pf * (X[-2, 1:] - X[-2, 0:-1])))
            term3 = -casadi.sum2(
                ds_vect.T
                * casadi.sum1((pa * uncertain_vals[0:size] - pf * kappa) * X[0:size, :])
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
        print("solving the Outer Optimization problem")
        start_time = time.time()
        sol = opti.solve()
        print(f"Done; time taken {time.time()-start_time} seconds...")

        # if _DEBUG:
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

        if scale == 1.0:
            print("Scaling the uncertainty values")
            theta_vals = uncertain_vals[:size]
            gamma_vals = uncertain_vals[size:]
            theta_vals_adj = theta_vals / thetaSD
            gamma_vals_adj = gamma_vals / gammaSD
            uncertain_vals_adj = np.concatenate((theta_vals_adj, gamma_vals_adj))
            site_stdev_adj = np.concatenate(
                (np.ones(thetaSD.shape), np.ones(gammaSD.shape))
            )
            site_covariances_adj = np.diag(np.power(site_stdev_adj, 2))
            site_precisions_adj = np.linalg.inv(site_covariances_adj)
            uncertain_vals_mean_adj = uncertain_vals_adj.copy()

            print("Uncertain_vals used in the sampling: ", uncertain_vals_adj)
            print("site_stdev_adj: ", site_stdev_adj)

            log_density = lambda uncertain_vals_adj: log_density_function(
                uncertain_vals=uncertain_vals_adj,
                uncertain_vals_mean=uncertain_vals_mean_adj,
                theta_vals=theta_vals_adj,
                thetaSD=thetaSD,
                gammaSD=gammaSD,
                scale=1.0,
                site_precisions=site_precisions_adj,
                alpha=alpha,
                N=N,
                # sol=sol,
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
                two_param_uncertainty=two_param_uncertainty,
            )
        else:
            print("Not scaling the uncertainty values")
            ## Start Sampling
            # Update signature of log density evaluator
            log_density = lambda uncertain_vals: log_density_function(
                uncertain_vals=uncertain_vals,
                uncertain_vals_mean=uncertain_vals_mean,
                theta_vals=theta_vals,
                thetaSD=None,
                gammaSD=None,
                scale=0.0,
                site_precisions=site_precisions,
                alpha=alpha,
                N=N,
                # sol=sol,
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
                two_param_uncertainty=two_param_uncertainty,
            )

        # Create MCMC sampler & sample, then calculate diagnostics
        if mode in [1.0, 3.0]:
            constraint_test_mode = lambda x: True if np.all(x >= 0) else False
        elif mode == 2.0:
            constraint_test_mode = lambda x: True if np.max(x >= 0) else False
        else:
            raise ValueError("Invalid mode")

        sampler = create_hmc_sampler(
            size=tot_size,
            log_density=log_density,
            #
            burn_in=100,
            mix_in=mix_in,
            symplectic_integrator="verlet",
            symplectic_integrator_stepsize=stepsize,
            symplectic_integrator_num_steps=symplectic_integrator_num_steps,
            constraint_test=constraint_test_mode,
            mass_matrix_gamma=mass_matrix_gamma,
            mass_matrix_theta=mass_matrix_theta,
            mass_matrix_gamma_scale=mass_matrix_gamma_scale,
            mass_matrix_theta_scale=mass_matrix_theta_scale,
        )

        if scale == 1.0:
            # Update to get the mode as well as the sample
            sampling_results = sampler.start_MCMC_sampling(
                sample_size=sample_size,
                initial_state=uncertain_vals_adj,
                verbose=True,
            )
        else:
            sampling_results = sampler.start_MCMC_sampling(
                sample_size=sample_size,
                initial_state=uncertain_vals,
                verbose=True,
            )

        count_positive = sampling_results["count_positive"]
        count_Ratio = count_positive / 2100
        print("positive_ratio", count_Ratio)
        count_Ratio_tracker.append(count_Ratio)

        site_count_positive = sampling_results["site_count_positive"]
        site_count_Ratio = site_count_positive / 2100
        print("site positive ratio", site_count_Ratio)
        site_count_Ratio_tracker.append(site_count_Ratio)

        uncertainty_post_samples = np.asarray(sampling_results["collected_ensemble"])
        uncertainty_map_estimate = sampling_results["map_estimate"]

        # Update ensemble/tracker
        collected_ensembles.update({cntr: uncertainty_post_samples.copy()})

        if mode_as_solution:
            print("uncertain values from last iteration: ", uncertain_vals_old)
            print("uncertain values from this iteration: ", uncertainty_map_estimate)
            if scale == 1.0:
                uncertainty_map_estimate[:size] *= thetaSD
                uncertainty_map_estimate[size:] *= gammaSD
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
            if scale == 1.0:
                uncertainty_post_samples[:, :size] *= thetaSD
                uncertainty_post_samples[:, size:] *= gammaSD
                print(
                    "scale back uncertain values from this iteration: ",
                    np.mean(uncertainty_post_samples, axis=0),
                )

            if mode == 2.0:  # drop negative values
                means = np.mean(uncertainty_post_samples, axis=0)
                stds = np.std(uncertainty_post_samples, axis=0)
                print(uncertainty_post_samples.shape)
                lower_bound = np.zeros(means.shape)
                upper_bound = np.inf * np.ones(means.shape)
                a, b = (lower_bound - means) / stds, (upper_bound - means) / stds
                trunc_normal_means = []
                trunc_normal_stds = []
                for i in range(uncertainty_post_samples.shape[1]):
                    distribution = truncnorm(a[i], b[i], loc=means[i], scale=stds[i])
                    trunc_mean = distribution.mean()
                    trunc_normal_means.append(trunc_mean)
                    trunc_std = distribution.std()
                    trunc_normal_stds.append(trunc_std)

                trunc_normal_means = np.array(trunc_normal_means)
                uncertain_vals = (
                    weight * trunc_normal_means + (1 - weight) * uncertain_vals_old
                )
                trunc_normal_stds = np.array(trunc_normal_stds)

                print("trunc_normal_means", trunc_normal_means)
                print("trunc_normal_stds", trunc_normal_stds)

                # mask = uncertainty_post_samples < 0
                # masked_samples = ma.masked_array(uncertainty_post_samples, mask)
                # mean_non_negative = np.mean(masked_samples, axis=0).data
                # print("mean_non_negative",mean_non_negative)
                # print("mean_with_negative",np.mean(uncertainty_post_samples, axis=0 ))
                # uncertain_vals = weight * np.mean(masked_samples, axis=0).data + (1-weight) * uncertain_vals_old
            else:
                uncertain_vals = (
                    weight * np.mean(uncertainty_post_samples, axis=0)
                    + (1 - weight) * uncertain_vals_old
                )

            print("updated uncertain values: ", uncertain_vals)

        uncertain_vals_tracker.append(uncertain_vals.copy())

        if mode == 2.0:
            uncertain_post_SD = trunc_normal_stds.copy()
        else:
            uncertain_post_SD = np.std(uncertainty_post_samples, axis=0)
        thetaSD = uncertain_post_SD[:size]
        gammaSD = uncertain_post_SD[size:]
        print(
            "uncertain value standard deviation from this iteration: ",
            uncertain_post_SD,
        )
        uncertain_SD_tracker.append(uncertain_post_SD.copy())

        mass_matrix = mass_matrix_weight * 1 / (uncertain_post_SD**2) + (
            1 - mass_matrix_weight
        ) * np.concatenate((mass_matrix_theta, mass_matrix_gamma))
        print("updated mass matrix:", mass_matrix)
        mass_matrix_theta = mass_matrix[0:size]
        mass_matrix_gamma = mass_matrix[size:]
        mass_matrix_tracker.append(mass_matrix.copy())

        # Evaluate error for convergence check
        # The percentage difference are changed to absolute difference
        abs_error = np.max(np.abs(uncertain_vals_old - uncertain_vals))
        percentage_error = np.max(
            np.abs(uncertain_vals_old - uncertain_vals) / uncertain_vals_old
        )
        log_diff_error = np.max(
            np.abs(np.log(uncertain_vals_old) - np.log(uncertain_vals))
        )

        abs_error_tracker.append(abs_error)
        percentage_error_tracker.append(percentage_error)
        log_diff_error_tracker.append(log_diff_error)

        print(
            decorate_text(
                f"Iteration [{cntr+1:4d}]: Absolte Error = {abs_error}, Percentage Error = {percentage_error}, Log Difference Error = {log_diff_error}"
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
                "count_Ratio_tracker": count_Ratio_tracker,
                "site_count_Ratio_tracker": site_count_Ratio_tracker,
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

    print("Terminated. Sampling the final distribution")
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


def solve_with_gams(
    # Configurations/Settings
    site_num=25,  # Number of sites(10, 25, 100, 1000)
    norm_fac=1e9,
    delta_t=0.02,
    alpha=0.045007414,
    kappa=2.094215255,
    pf=20.76,
    pa=44.75,
    xi=0.01,
    zeta=1.66e-4 * 1e9,  # zeta := 1.66e-4*norm_fac  #
    #
    max_iter=20000,
    tol=0.001,
    T=200,
    N=200,
    #
    sample_size=1000,  # simulations before convergence (to evaluate the mean)
    mode_as_solution=False,  # If true, use the mode (point of high probability) as solution for gamma
    final_sample_size=25_000,  # number of samples to collect after convergence
    two_param_uncertainty=True,
    weight=0.25,  # <-- Not sure how this linear combination weighting helps!
    output_dir="GAMS_Results",
    mix_in=2,
    mass_matrix_theta_scale=1.0,
    mass_matrix_gamma_scale=1.0,
    symplectic_integrator_num_steps=10,
    mass_matrix_weight=0.1,
    stepsize=0.1,
):
    """
    Main function to solve the bilievel optimization problem using gams for the outer optimization problem.

    :param float tol: convergence tolerance
    :param T:
    :param N:
    """
    # if gams is None:
    #     print("Failed to import GAMS. This function requires GAMSto be installed!")
    #     raise ImportError

    # Create the output directory
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # Load sites' data
    (
        zbar_2017,
        gamma,
        gammaSD,
        z_2017,
        forestArea_2017_ha,
        theta,
        thetaSD,
    ) = load_site_data(
        site_num,
        norm_fac=norm_fac,
    )

    # gammaSD_adj = np.ones(gamma.shape)
    # gamma_adj=gammaSD_adj/gammaSD
    # gamma=gamma*gamma_adj
    # gammaSD=gammaSD_adj.copy()
    # thetaSD_adj = np.ones(theta.shape)
    # theta_adj=thetaSD_adj/thetaSD
    # theta=theta*theta_adj
    # thetaSD=thetaSD_adj.copy()
    print("theta", theta, "tehtaSD", thetaSD)
    print("gamma", gamma, "gammaSD", gammaSD)

    df_z_2017 = pd.DataFrame(z_2017)
    saveto = os.path.join(output_dir, "z0Data.csv")
    df_z_2017.to_csv(saveto)
    df_zbar_2017 = pd.DataFrame(zbar_2017)
    saveto = os.path.join(output_dir, "zbarData.csv")
    df_zbar_2017.to_csv(saveto)

    # Evaluate Gamma values ()
    gamma_1_vals = gamma - gammaSD
    gamma_2_vals = gamma + gammaSD
    gamma_vals = gamma
    size = gamma.size

    # Theta Values
    theta_vals = theta
    tot_size = gamma.size + theta.size

    # Retrieve z data for selected site(s)
    site_z_vals = z_2017

    if not two_param_uncertainty:
        # One parameter (gamma) uncertainty

        # Evaluate mean and covariances from site data
        site_stdev = gammaSD
        site_covariances = np.diag(np.power(site_stdev, 2))
        site_precisions = np.linalg.inv(site_covariances)
        site_mean = gamma_1_vals / 2 + gamma_2_vals / 2

        # Initialize Gamma Values
        uncertain_vals = gamma.copy()
        uncertain_vals_mean = gamma.copy()
        uncertain_vals_old = gamma.copy()

    else:
        # Two parameter uncertainty (both theta and gamma)

        vals = np.concatenate((theta_vals, gamma_vals))
        # Evaluate mean and covariances from site data
        site_stdev = np.concatenate((thetaSD, gammaSD))
        site_covariances = np.diag(np.power(site_stdev, 2))
        site_precisions = np.linalg.inv(site_covariances)
        site_mean = vals

        # Initialize Gamma Values
        uncertain_vals = vals.copy()
        uncertain_vals_mean = vals.copy()
        uncertain_vals_old = vals.copy()

    # Householder to track sampled gamma values
    # uncertain_vals_tracker       = np.empty((uncertain_vals.size, sample_size+1))
    # uncertain_vals_tracker[:, 0] = uncertain_vals.copy()
    uncertain_vals_tracker = [uncertain_vals.copy()]
    uncertain_SD_tracker = [np.concatenate((thetaSD, gammaSD)).copy()]

    mass_matrix_theta = thetaSD**2
    mass_matrix_gamma = gammaSD**2
    mass_matrix = np.concatenate((mass_matrix_theta, mass_matrix_gamma))
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
        mass_matrix_theta_scale=mass_matrix_theta_scale,
        mass_matrix_gamma_scale=mass_matrix_gamma_scale,
        symplectic_integrator_num_steps=symplectic_integrator_num_steps,
        two_param_uncertainty=two_param_uncertainty,
        gammaSD=gammaSD,
        thetaSD=thetaSD,
        weight=weight,
        mass_matrix_weight=mass_matrix_weight,
        output_dir=output_dir,
        stepsize=stepsize,
    )

    # Initialize error & iteration counter
    abs_error = np.infty
    percentage_error = np.infty
    log_diff_error = np.infty
    cntr = 0

    # Loop until convergence
    while cntr < max_iter and percentage_error > tol:
        print(decorate_text(f"Optimization Iteration[{cntr+1}/{max_iter}]"))

        if not two_param_uncertainty:
            # One parameter (gamma) uncertainty

            # Update x0
            x0_vals = uncertain_vals * forestArea_2017_ha

            x0data = pd.DataFrame(x0_vals)
            saveto = os.path.join(output_dir, "X0Data.csv")
            x0data.to_csv(saveto)

            gammadata = pd.DataFrame(uncertain_vals)
            saveto = os.path.join(output_dir, "GammaData.csv")
            gammadata.to_csv(saveto)

            # Create Gams Workspace
            username = getpass.getuser()
            if username == "hqin":
                ws = GamsWorkspace(
                    system_directory="/home/hqin/gams/gams43.4_linux_x64_64_sfx/",
                    working_directory=output_dir,
                )
            elif username == "pengyu":
                ws = GamsWorkspace(
                    system_directory="/home/pengyu/gams43.4_linux_x64_64_sfx",
                    working_directory=output_dir,
                )

            # TODO: I am not sure where these GAMS model files are generated!
            # We may need some gams_model_dir to put these files in for other sites as well!
            # print(f"amazon_{size}sites.gms")
            t1 = ws.add_job_from_file(f"amazon_{size}sites.gms")
            t1.run()

            readfrom = os.path.join(output_dir, "amazon_data_u.dat")
            dfu = pd.read_csv(readfrom, delimiter="\t")
            # Process the data using the pandas DataFrame
            dfu = dfu.drop("T/R ", axis=1)
            sol_val_Up = dfu.to_numpy()

            readfrom = os.path.join(output_dir, "amazon_data_w.dat")
            dfw = pd.read_csv(readfrom, delimiter="\t")
            # Process the data using the pandas DataFrame
            dfw = dfw.drop("T   ", axis=1)
            dfw_np = dfw.to_numpy()

            readfrom = os.path.join(output_dir, "amazon_data_x.dat")
            dfx = pd.read_csv(readfrom, delimiter="\t")
            # Process the data using the pandas DataFrame
            dfx = dfx.drop("T   ", axis=1)
            dfx_np = dfx.to_numpy()

            readfrom = os.path.join(output_dir, "amazon_data_z.dat")
            dfz = pd.read_csv(readfrom, delimiter="\t")
            # Process the data using the pandas DataFrame
            dfz = dfz.drop("T/R ", axis=1)
            dfz_np = dfz.to_numpy()

            sol_val_Ua = dfw_np**2
            sol_val_X = np.concatenate((dfz_np.T, dfx_np.T))

        else:
            # Two parameter uncertainty (both theta and gamma)

            # Update x0
            x0_vals = uncertain_vals[size:] * forestArea_2017_ha

            x0data = pd.DataFrame(x0_vals)
            saveto = os.path.join(output_dir, "X0Data.csv")
            x0data.to_csv(saveto)

            gammadata = pd.DataFrame(uncertain_vals[size:])
            saveto = os.path.join(output_dir, "GammaData.csv")
            gammadata.to_csv(saveto)

            thetadata = pd.DataFrame(uncertain_vals[0:size])
            saveto = os.path.join(output_dir, "ThetaData.csv")
            thetadata.to_csv(saveto)

            # Create Gams Workspace
            username = getpass.getuser()
            if username == "hqin":
                ws = GamsWorkspace(
                    system_directory="/home/hqin/gams/gams43.4_linux_x64_64_sfx/",
                    working_directory=output_dir,
                )

            elif username == "pengyu":
                ws = GamsWorkspace(
                    system_directory="/home/pengyu/gams43.4_linux_x64_64_sfx",
                    working_directory=output_dir,
                )
            print("GAMS workspace created: " + username)

            gams_file = f"amazon_{size}sites.gms"
            shutil.copy(gams_file, output_dir)
            t1 = ws.add_job_from_file(
                os.path.join(output_dir, os.path.basename(gams_file))
            )
            t1.run()

            readfrom = os.path.join(output_dir, "amazon_data_u.dat")
            dfu = pd.read_csv(readfrom, delimiter="\t")

            # Process the data using the pandas DataFrame
            dfu = dfu.drop("T/R ", axis=1)
            sol_val_Up = dfu.to_numpy()

            readfrom = os.path.join(output_dir, "amazon_data_w.dat")
            dfw = pd.read_csv(readfrom, delimiter="\t")
            # Process the data using the pandas DataFrame
            dfw = dfw.drop("T   ", axis=1)
            dfw_np = dfw.to_numpy()

            readfrom = os.path.join(output_dir, "amazon_data_x.dat")
            dfx = pd.read_csv(readfrom, delimiter="\t")
            # Process the data using the pandas DataFrame
            dfx = dfx.drop("T   ", axis=1)
            dfx_np = dfx.to_numpy()

            readfrom = os.path.join(output_dir, "amazon_data_z.dat")
            dfz = pd.read_csv(readfrom, delimiter="\t")
            # Process the data using the pandas DataFrame
            dfz = dfz.drop("T/R ", axis=1)
            dfz_np = dfz.to_numpy()

            sol_val_Ua = dfw_np**2
            sol_val_X = np.concatenate((dfz_np.T, dfx_np.T))

        print("sol.value(X)", sol_val_X)
        print("sol.value(Ua)", sol_val_Ua)
        print("sol.value(Up)", sol_val_Up)

        sol_val_X_tracker.append(sol_val_X)
        sol_val_Ua_tracker.append(sol_val_Ua)
        sol_val_Up_tracker.append(sol_val_Up)
        # sol_val_Um_tracker.append(sol_val_Um)
        # sol_val_Z_tracker.append(sol_val_Z)

        ## Start Sampling
        # Update signature of log density evaluator
        log_density = lambda uncertain_vals: log_density_function(
            uncertain_vals=uncertain_vals,
            uncertain_vals_mean=uncertain_vals_mean,
            theta_vals=theta_vals,
            site_precisions=site_precisions,
            alpha=alpha,
            N=N,
            # sol=sol,
            sol_val_X=sol_val_X,
            sol_val_Ua=sol_val_Ua[0:-1].T,
            sol_val_Up=sol_val_Up[0:-1].T,
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
            two_param_uncertainty=two_param_uncertainty,
        )

        # Create MCMC sampler & sample, then calculate diagnostics
        sampler = create_hmc_sampler(
            size=tot_size,
            log_density=log_density,
            #
            burn_in=100,
            mix_in=mix_in,
            symplectic_integrator="verlet",
            symplectic_integrator_stepsize=stepsize,
            symplectic_integrator_num_steps=symplectic_integrator_num_steps,
            constraint_test=lambda x: True if np.all(x >= 0) else False,
            mass_matrix_gamma=mass_matrix_gamma,
            mass_matrix_theta=mass_matrix_theta,
            mass_matrix_gamma_scale=mass_matrix_gamma_scale,
            mass_matrix_theta_scale=mass_matrix_theta_scale,
        )

        # Update to get the mode as well as the sample
        sampling_results = sampler.start_MCMC_sampling(
            sample_size=sample_size,
            initial_state=uncertain_vals,
            verbose=True,
        )
        uncertainty_post_samples = np.asarray(sampling_results["collected_ensemble"])
        uncertainty_map_estimate = sampling_results["map_estimate"]

        # Update ensemble/tracker
        collected_ensembles.update({cntr: uncertainty_post_samples.copy()})

        # Update gamma value
        if mode_as_solution:
            uncertain_vals = (
                weight * uncertainty_map_estimate + (1 - weight) * uncertain_vals_old
            )

        else:
            uncertain_vals = (
                weight * np.mean(uncertainty_post_samples, axis=0)
                + (1 - weight) * uncertain_vals_old
            )
        uncertain_vals_tracker.append(uncertain_vals.copy())

        uncertain_post_SD = np.std(uncertainty_post_samples, axis=0)
        thetaSD = uncertain_post_SD[0:size]
        gammaSD = uncertain_post_SD[size:]
        uncertain_SD_tracker.append(uncertain_post_SD.copy())

        mass_matrix = mass_matrix_weight * uncertain_post_SD**2 + (
            1 - mass_matrix_weight
        ) * np.concatenate((mass_matrix_theta, mass_matrix_gamma))
        mass_matrix_theta = mass_matrix[0:size]
        mass_matrix_gamma = mass_matrix[size:]
        mass_matrix_tracker.append(mass_matrix.copy())

        # Evaluate error for convergence check
        # The percentage difference are changed to absolute difference
        abs_error = np.max(np.abs(uncertain_vals_old - uncertain_vals))
        percentage_error = np.max(
            np.abs(uncertain_vals_old - uncertain_vals) / uncertain_vals_old
        )
        log_diff_error = np.max(
            np.abs(np.log(uncertain_vals_old) - np.log(uncertain_vals))
        )

        abs_error_tracker.append(abs_error)
        percentage_error_tracker.append(percentage_error)
        log_diff_error_tracker.append(log_diff_error)
        print(
            decorate_text(
                f"Iteration [{cntr+1:4d}]: Absolte Error = {abs_error}, Percentage Error = {percentage_error}, Log Difference Error = {log_diff_error}"
            )
        )

        # Exchange gamma values (for future weighting/update & error evaluation)
        uncertain_vals_old = uncertain_vals

        # Increase the counter
        cntr += 1

        results.update(
            {
                "cntr": cntr,
                "abs_error_tracker": np.asarray(abs_error_tracker),
                "percentage_error_tracker": np.asarray(percentage_error_tracker),
                "log_diff_error_tracker": np.asarray(log_diff_error_tracker),
                "uncertain_vals_tracker": np.asarray(uncertain_vals_tracker),
                "uncertain_SD_tracker": np.asarray(uncertain_SD_tracker),
                "collected_ensembles": collected_ensembles,
                "sol_val_X_tracker": sol_val_X_tracker,
                "sol_val_Ua_tracker": sol_val_Ua_tracker,
                "sol_val_Up_tracker": sol_val_Up_tracker,
                "sol_val_Um_tracker": sol_val_Um_tracker,
                "sol_val_Z_tracker": sol_val_Z_tracker,
            }
        )

        saveto = os.path.join(output_dir, "results.pcl")
        pickle.dump(results, open(saveto, "wb"))
        print(f"Results saved to {saveto}")

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

    print("Terminated. Sampling the final distribution")
    # Sample (densly) the final distribution
    final_sample = sampler.sample(
        sample_size=final_sample_size,
        initial_state=uncertain_vals,
        verbose=True,
    )
    final_sample = np.asarray(final_sample)
    results.update({"final_sample": final_sample})
    saveto = os.path.join(output_dir, "results.pcl")
    pickle.dump(results, open(saveto, "wb"))
    print(f"Results saved to {saveto}")

    return results
