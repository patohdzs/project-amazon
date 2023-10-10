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

import pysrc.plots as plots
from pysrc.services.data_service import load_site_data
from pysrc.services.file_service import (
    get_path,
    logs_dir_path,
    output_dir_path,
    plots_dir_path,
)
from pysrc.solvers import gamma_fitted, theta_fitted

sys.path.append(os.path.abspath("src"))
sns.set(font_scale=1.2)

# Read arguments from stdin
parser = argparse.ArgumentParser(description="parameter settings")

parser.add_argument("--model", type=str, default="calibrated_coef_priors")
parser.add_argument("--xi", type=float, default=1.0)
parser.add_argument("--pf", type=float, default=25)
parser.add_argument("--pa", type=float, default=44.75)
parser.add_argument("--weight", type=float, default=0.25)
parser.add_argument("--sitenum", type=int, default=10)
parser.add_argument("--timehzn", type=int, default=200)

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
print("Transforming prior samples...")
beta_theta_prior_samples = pd.read_csv(
    get_path("data", "hmc", "theta_coe_ori.csv")
).to_numpy()[:5000, :]

beta_gamma_prior_samples = pd.read_csv(
    get_path("data", "hmc", "gamma_coe_ori.csv")
).to_numpy()

theta_prior_samples = np.array(
    [theta_fitted(c, site_theta_2017_df) for c in beta_theta_prior_samples]
)
gamma_prior_samples = np.array(
    [gamma_fitted(c, site_gamma_2017_df) for c in beta_gamma_prior_samples]
)

prior_samples = np.concatenate((theta_prior_samples, gamma_prior_samples), axis=1)
print("Done!")

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
