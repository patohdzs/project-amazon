import getpass
import os
import pickle
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mcmc.hmc import create_hmc_sampler
from services.data_service import load_site_data
from solvers import log_density_function
from utils.text import decorate_text

# Check if gams is available; delay exception raise to the call
try:
    from gams import GamsWorkspace

    if GamsWorkspace.api_major_rel_number < 42:  # old API structure
        from gams import *  # noqa: F403
    else:  # new API structure
        from gams.control import *  # noqa: F403
except ImportError:
    gams = None


_GAMS_SYSTEM_LOADER = os.path.join(
    os.path.dirname(__file__), "_gams_system_directory.dat"
)


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
        prompt = "\n**\nPlease input FULL path to GAMS system directory/resources.\n"
        prompt += """For example:
                    '/Library/Frameworks/GAMS.framework/Versions/43/Resources'\n"""
        gams_sys_dir = input(prompt).strip(""" \n" '  """)
        # Write it to file and recurse
        with open(filepath, "w") as f_id:
            f_id.write(gams_sys_dir)

        # Recurse to validate the path
        return get_gams_system_directory(
            filepath=filepath,
        )

    return gams_sys_dir


def solve_with_gams(
    # Configurations/Settings
    site_num=25,  # Number of sites(10, 25, 100, 1000)
    norm_fac=1e12,
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
    mode_as_solution=False,  # If true, use the mode as solution for gamma
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
    Main function to solve the bilievel optimization problem using
    gams for the outer optimization problem.

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

    if not two_param_uncertainty:
        # One parameter (gamma) uncertainty

        # Evaluate mean and covariances from site data
        site_stdev = gammaSD
        site_covariances = np.diag(np.power(site_stdev, 2))
        site_precisions = np.linalg.inv(site_covariances)
        gamma_1_vals / 2 + gamma_2_vals / 2

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

        # Initialize Gamma Values
        uncertain_vals = vals.copy()
        uncertain_vals_mean = vals.copy()
        uncertain_vals_old = vals.copy()

    # Householder to track sampled gamma values
    # uncertain_vals_tracker = np.empty((uncertain_vals.size, sample_size+1))
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
            # We may need gams_model_dir to put these files in for other sites as well!
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
        def log_density(uncertain_vals):
            return log_density_function(
                uncertain_vals=uncertain_vals,
                uncertain_vals_mean=uncertain_vals_mean,
                theta_vals=theta_vals,
                site_precisions=site_precisions,
                alpha=alpha,
                N=N,
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
            mass_matrix=mass_matrix,
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
                f"""
                Iteration [{cntr+1:4d}]: Absolte Error = {abs_error},
                Percentage Error = {percentage_error},
                Log Difference Error = {log_diff_error}
                """
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
