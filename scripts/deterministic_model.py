import argparse

from pysrc.optimization import gurobi
from pysrc.services.data_service import load_site_data

# Read arguments from stdin
parser = argparse.ArgumentParser(description="parameter settings")

parser.add_argument("--model", type=str, default="deterministic")
parser.add_argument("--xi", type=float, default=1.0)
parser.add_argument("--pf", type=float, default=25)
parser.add_argument("--pa", type=float, default=44.75)
parser.add_argument("--weight", type=float, default=0.25)
parser.add_argument("--sitenum", type=int, default=10)
parser.add_argument("--timehzn", type=int, default=200)

# Parse arguments
args = parser.parse_args()

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

(
    sol_val_X,
    sol_val_Up,
    sol_val_Um,
    sol_val_Z,
    sol_val_Ua,
) = gurobi.solve_planner_problem(
    T=args.timehzn,
    theta=theta,
    gamma=gamma,
    x0=x0_vals,
    zbar=zbar_2017,
    z0=z_2017,
    pe=args.pf,
    pa=args.pa,
)
