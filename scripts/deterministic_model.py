import argparse
import pickle

from pysrc.optimization import gams, gurobi
from pysrc.services.data_service import load_site_data
from pysrc.services.file_service import logs_dir_path, output_dir_path, plots_dir_path

# Read arguments from stdin
parser = argparse.ArgumentParser(description="parameter settings")

parser.add_argument("--model", type=str, default="deterministic")
parser.add_argument("--opt", type=str, default="gurobi")
parser.add_argument("--pe", type=float, default=25)
parser.add_argument("--pa", type=float, default=41.11)
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
    gamma,
    z_2017,
    forest_area_2017,
    theta,
    _,
    _,
    _,
    _,
) = load_site_data(args.sitenum)

# Computing carbon absorbed in start period
x0_vals = gamma * forest_area_2017

# Choose optimizer
if args.opt == "gurobi":
    solve_planner_problem = gurobi.solve_planner_problem

elif args.opt == "gams":
    solve_planner_problem = gams.solve_planner_problem

else:
    raise ValueError("Optimizer must be one of ['gurobi', 'gams']")

results = solve_planner_problem(
    T=args.timehzn,
    theta=theta,
    gamma=gamma,
    x0=x0_vals,
    zbar=zbar_2017,
    z0=z_2017,
    pe=args.pe,
    pa=args.pa,
)

# Save results
outfile_path = output_dir / "results.pcl"
with open(outfile_path, "wb") as outfile:
    pickle.dump(results, outfile)
    print(f"Results saved to {outfile_path}")
