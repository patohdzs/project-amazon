import argparse
import os
import pickle

from pysrc import plots
from pysrc.sampling import baseline
from pysrc.services.file_service import plots_dir_path

# Read arguments from stdin
parser = argparse.ArgumentParser(description="parameter settings")

parser.add_argument("--model", type=str, default="full_model")
parser.add_argument("--sitenum", type=int, default=10)

# Parse arguments
args = parser.parse_args()

# Create plots directories
plots_dir = plots_dir_path(**vars(args))

# Sample baseline
fit = baseline.sample(args.model, 4000, args.sitenum)

theta = fit.stan_variable("theta").T
gamma = fit.stan_variable("gamma").T
# Directory to save results
results_dir = "results"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

filename = "prior.pcl"
file_path = os.path.join(results_dir, filename)

data = {"theta": theta, "gamma": gamma}

with open(file_path, "wb") as f:
    pickle.dump(data, f)

print(f"Results (theta and gamma) saved to {file_path}")

# Plotting theta and gamma
plots.basline_density(
    fit["theta"].T, fit["gamma"].T, num_sites=args.sitenum, plots_dir=plots_dir
)
