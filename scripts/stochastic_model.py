import argparse
import pickle

import numpy as np

from pysrc.optimization.mpc import price_path_probs, price_paths, solve_planner_problem
from pysrc.sampling import baseline
from pysrc.sampling.hmm_estimation import estimate_price_model
from pysrc.sampling.mpc import sample_price_path
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
price_states, M = estimate_price_model()

# Expand initial conditions
z0 = np.tile(np.reshape(z_2017, (-1, 1)), (1, 2**args.tau))
x0 = np.tile(np.reshape(x_2017, (-1, 1)), (1, 2**args.tau))

# Compute sample price path
sample_pa_path = sample_price_path(args.timehzn, M, price_states)
full_trajectory = []
for pa in sample_pa_path[:3]:
    # Determine if currenctly on high price
    start_high = pa == price_states[1]

    # Compute possible price paths going forward
    pa_paths = price_paths(args.timehzn, args.tau, price_states)

    # Compute likelihood of each path
    pa_path_probs = price_path_probs(M, args.tau, pa_paths)

    # Solve planner problem with stochastic horizon
    trajectories = solve_planner_problem(
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

    traj = {
        "Z": trajectories["Z"][1, :, 1],
        "X": trajectories["X"][1, :, 1],
        "U": trajectories["U"][1, :, 1],
        "V": trajectories["V"][1, :, 1],
        "w": trajectories["w"][1, 1],
    }

    # Save full trajectory
    full_trajectory.append(traj)

    # Re-set initial conditions
    z0 = np.tile(np.reshape(traj["Z"], (-1, 1)), (1, 2**args.tau))
    x0 = np.tile(np.reshape(traj["X"], (-1, 1)), (1, 2**args.tau))


# Save results
outfile_path = output_dir / "results.pcl"
with open(outfile_path, "wb") as outfile:
    pickle.dump(full_trajectory, outfile)
    print(f"Results saved to {outfile_path}")
