import argparse
import os
import pickle
import sys

import numpy as np
import seaborn as sns

from pysrc import plots
from pysrc.sampling import baseline
from pysrc.services.data_service import load_site_data
from pysrc.services.file_service import (
    logs_dir_path,
    output_dir_path,
    plots_dir_path,
)

sys.path.append(os.path.abspath("src"))
sns.set(font_scale=1.2)

# Read arguments from stdin
parser = argparse.ArgumentParser(description="parameter settings")

parser.add_argument("--model", type=str, default="full_model_v3")
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


with open(output_dir / "results.pcl", "rb") as f:
    # Load the data from the file
    results = pickle.load(f)

# Load site data
(
    zbar_2017,
    gamma,
    z_2017,
    forestArea_2017_ha,
    theta,
    site_theta_2017_df,
    site_gamma_2017_df,
    municipal_theta_df,
    municipal_gamma_df,
) = load_site_data(args.sitenum)


# Load coef baseline samples
fit = baseline.sample(model_name=args.model, num_samples=10000, num_sites=args.sitenum)
base_samples = np.concatenate(
    (fit.stan_variable("theta"), fit.stan_variable("gamma")), axis=1
)

try:
    adj_samples = results["final_sample"]
except KeyError:
    print(
        """
        Algorithm did not finish converging.
        Plotting adjusted samples from last iteration...
        """
    )
    adj_samples = results["collected_ensembles"][results["cntr"] - 1]

# Plot overlapped densities
plots.density_overlap(base_samples, adj_samples, plots_dir, args.sitenum)

# Plot absolute and percentage error
plots.traceplot_abs_error(results, plots_dir)
plots.traceplot_pct_error(results, plots_dir)

# Plot sampling time
plots.traceplot_sampling_time(results, plots_dir)

# Plot theta and gamma convergence
plots.traceplot_params_pct_error(results, plots_dir)
plots.traceplot_gammas(results, plots_dir)
plots.traceplot_thetas(results, plots_dir)

# Plot agg %change in z
plots.agg_Z_trajectory(
    z_2017=z_2017,
    zbar_2017=zbar_2017,
    results=results,
    num_sites=args.sitenum,
    plots_dir=plots_dir,
)

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
