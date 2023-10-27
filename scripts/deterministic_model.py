import argparse

from pysrc.optimization import gurobi
from pysrc.services.data_service import load_site_data

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

# Load site data
(
    zbar_2017,
    gamma,
    z_2017,
    forestArea_2017_ha,
    theta,
    _,
    _,
    _,
    _,
) = load_site_data(args.sitenum)

# Computing carbon absorbed in start period
x0_vals = gamma * forestArea_2017_ha

gurobi.solve_outer_optimization_problem(
    T=args.timehzn,
    theta=theta,
    gamma=gamma,
    x0_vals=x0_vals,
    zbar_2017=zbar_2017,
    z_2017=z_2017,
)
