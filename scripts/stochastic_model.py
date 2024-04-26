import argparse
import pickle

import numpy as np

from pysrc.optimization.mpc import price_path_probs, price_paths, solve_planner_problem
from pysrc.sampling import baseline
from pysrc.sampling.hmm_estimation import estimate_price_model
from pysrc.services.data_service import load_site_data
from pysrc.services.file_service import logs_dir_path, output_dir_path, plots_dir_path

# Read arguments from stdin
parser = argparse.ArgumentParser(description="parameter settings")

parser.add_argument("--model", type=str, default="deterministic")
parser.add_argument("--opt", type=str, default="gurobi")
parser.add_argument("--pe", type=float, default=25)
parser.add_argument("--tau", type=float, default=3)
parser.add_argument("--sitenum", type=int, default=78)
parser.add_argument("--timehzn", type=int, default=200)

# Parse arguments
args = parser.parse_args()

# Create output and plots directories
output_dir = output_dir_path(**vars(args))
plots_dir = plots_dir_path(**vars(args))
logs_dir = logs_dir_path(**vars(args))

# Load site data
(
    zbar_2017,
    z_2017,
    forest_area_2017,
    _,
    _,
    _,
    _,
) = load_site_data(args.sitenum)

# Set productivity parameters using baseline mean
baseline_fit = baseline.sample(
    model_name="full_model",
    num_sites=args.sitenum,
    iter_sampling=10**4,
    chains=5,
    seed=1,
)

theta = baseline_fit.stan_variable("theta").mean(axis=0)
gamma = baseline_fit.stan_variable("gamma").mean(axis=0)

# Computing carbon absorbed in start period
x_2017 = gamma * forest_area_2017

# Estimate HMM model for prices
(
    aic,
    ll,
    bic,
    states,
    sigmas,
    P,
    sta_dist,
    sta_price,
    M,
) = estimate_price_model(
    [0.5, 0.5],
)


# Expand initial conditions
x0 = np.tile(np.reshape(x_2017, (-1, 1)), (1, 2**args.tau))
z0 = np.tile(np.reshape(z_2017, (-1, 1)), (1, 2**args.tau))
print(x0.shape)

# Compute possible price paths
pa_paths = price_paths(args.timehzn, args.tau, states)

# Compute likelihood of each path
pa_path_probs = price_path_probs(M, args.tau, pa_paths)

# Solve planner problem with stochastic horizon
results = solve_planner_problem(
    T=args.timehzn,
    theta=theta,
    gamma=gamma,
    x0=x0,
    z0=z0,
    zbar=zbar_2017,
    pa_paths=pa_paths,
    pa_path_probs=pa_path_probs,
    pe=args.pe,
)

print(results["Z"].shape)
print(results["Z"][1, :, :])
print(results["Z"][1, :, :].shape)


# Save results
outfile_path = output_dir / "results.pcl"
with open(outfile_path, "wb") as outfile:
    pickle.dump(results, outfile)
    print(f"Results saved to {outfile_path}")
