#!/usr/bin/env python

# Import Required Packages
# ========================
import argparse
import os
import pickle
import sys

import numpy as np
import pandas as pd
import seaborn as sns
from services.data_service import load_site_data
from services.file_service import logs_dir_path, output_dir_path, plots_dir_path
from solvers import gamma_fitted, theta_fitted

import plots

sys.path.append(os.path.abspath("src"))
sns.set(font_scale=1.2)

# Read arguments from stdin
parser = argparse.ArgumentParser(description="parameter settings")

parser.add_argument("--dataname", type=str, default="stan_lognorm")
parser.add_argument("--weight", type=float, default=0.25)
parser.add_argument("--xi", type=float, default=0.01)
parser.add_argument("--pf", type=float, default=20.76)
parser.add_argument("--pa", type=float, default=44.75)
parser.add_argument("--sitenum", type=int, default=10)
parser.add_argument("--time", type=int, default=200)
parser.add_argument("--mix_in", type=int, default=2)
parser.add_argument("--mass_matrix_theta_scale", type=float, default=1.0)
parser.add_argument("--mass_matrix_gamma_scale", type=float, default=1.0)
parser.add_argument("--mass_matrix_weight", type=float, default=0.1)
parser.add_argument("--symplectic_integrator_num_steps", type=int, default=2)
parser.add_argument("--stepsize", type=float, default=0.1)
parser.add_argument("--scale", type=float, default=1.0)
parser.add_argument("--mode", type=float, default=1.0)

# Parse arguments
args = parser.parse_args()

# Create output and plots directories
output_dir = output_dir_path(**vars(args))
plots_dir = plots_dir_path(**vars(args))
logs_dir = logs_dir_path(**vars(args))

# Load data
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
    gamma_coef_vcov,
    theta_coef_vcov,
    site_theta_2017_df,
    site_gamma_2017_df,
) = load_site_data(args.sitenum, norm_fac=1e11)


with open(output_dir / "results.pcl", "rb") as f:
    # Load the data from the file
    results = pickle.load(f)

# Prior samples
beta_theta_prior_samples = pd.read_csv("./data/hmc/theta_coe_ori.csv").to_numpy()
beta_gamma_prior_samples = pd.read_csv("./data/hmc/gamma_coe_ori.csv").to_numpy()

theta_prior_samples = np.array(
    [theta_fitted(c, site_theta_2017_df) for c in beta_theta_prior_samples]
)
gamma_prior_samples = np.array(
    [gamma_fitted(c, site_gamma_2017_df) for c in beta_gamma_prior_samples]
)

prior_samples = np.concatenate(
    (theta_prior_samples[:5000, :], gamma_prior_samples), axis=1
)


# Ploting historgam of prior samples
# plots.prior_density(samples=prior_samples, plots_dir=plots_dir, num_sites=10)


try:
    post_samples = results["final_sample"]
except KeyError:
    print(
        """
        Algorithm did not finish converging.
        Plotting posterior samples from last iteration...
        """
    )
    post_samples = results["collected_ensembles"][results["cntr"] - 1]

# plots.posterior_density(post_samples, plots_dir)


# Plot overlapped prior-posterior
plots.overlap_prior_posterior(
    prior_samples,
    post_samples,
    plots_dir,
)


# Plot absolute and percentage error
plots.traceplot_abs_error(results, plots_dir)
plots.traceplot_pct_error(results, plots_dir)


# Plot theta and gamma convergence
plots.traceplot_params_pct_error(results, plots_dir)
plots.traceplot_gammas(results, plots_dir)
plots.traceplot_thetas(results, plots_dir)

# Plot Z (for each site)
plots.Z_trajectory(results, plots_dir)

# Plot X (aggregate among sites)
plots.X_trajectory(results, plots_dir)

plots.delta_Z_trajectory(results, plots_dir)

# Plot adjustments
plots.adj_costs_trajectory(results, plots_dir)

# Plot trajectory of controls
plots.V_trajectory(results, plots_dir)
plots.U_trajectory(results, plots_dir)
