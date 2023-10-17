import argparse

from sampling.priors import sample_priors
from services.file_service import plots_dir_path

import plots

# Read arguments from stdin
parser = argparse.ArgumentParser(description="parameter settings")

parser.add_argument("--model", type=str, default="full_model_v3")
parser.add_argument("--sitenum", type=int, default=10)

# Parse arguments
args = parser.parse_args()

# Create plots directories
plots_dir = plots_dir_path(**vars(args))

# Sample prior
fit = sample_priors(args.model, 1000, args.sitenum)

# Plotting theta and gamma
plots.prior_density(
    fit["theta"].T, fit["gamma"].T, num_sites=args.sitenum, plots_dir=plots_dir
)
