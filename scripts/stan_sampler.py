#!/usr/bin/env python

# Import Required Packages
# ========================
import argparse

from pysrc.services.file_service import logs_dir_path, output_dir_path, plots_dir_path

# Import the solvers
from pysrc.solvers.stan import sample_with_stan

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

# Solve model with Casadi
stan_results = sample_with_stan(
    model_name=args.model + ".stan",
    output_dir=output_dir,
    xi=args.xi,
    pf=args.pf,
    pa=args.pa,
    weight=args.weight,
    site_num=args.sitenum,
    T=args.timehzn,
    max_iter=50,
    sample_size=1000,
    final_sample_size=5_000,
    num_chains=8,
)
