#!/usr/bin/env python

# Import Required Packages
# ========================
import argparse

from services.file_service import logs_dir_path, output_dir_path, plots_dir_path

# Import the solvers
from solvers.casadi import solve_with_casadi

# Read arguments from stdin
parser = argparse.ArgumentParser(description="parameter settings")

parser.add_argument("--dataname", type=str, default="tests")
parser.add_argument("--weight", type=float, default=0.25)
parser.add_argument("--xi", type=float, default=0.01)
parser.add_argument("--pf", type=float, default=20.76)
parser.add_argument("--pa", type=float, default=44.75)
parser.add_argument("--theta", type=float, default=1.0)
parser.add_argument("--gamma", type=float, default=1.0)
parser.add_argument("--sitenum", type=int, default=10)
parser.add_argument("--time", type=int, default=200)
parser.add_argument("--mix_in", type=int, default=2)
parser.add_argument("--mass_matrix_theta_scale", type=float, default=1.0)
parser.add_argument("--mass_matrix_gamma_scale", type=float, default=1.0)
parser.add_argument("--mass_matrix_weight", type=float, default=0.1)
parser.add_argument("--symplectic_integrator_num_steps", type=int, default=2)
parser.add_argument("--stepsize", type=float, default=0.1)
parser.add_argument("--scale", type=float, default=1.0)
parser.add_argument("--mode", type=float, default=1.0)

# Parse arguments
args = parser.parse_args()

# Create output and plots directories
output_dir = output_dir_path(**vars(args))
plots_dir = plots_dir_path(**vars(args))
logs_dir = logs_dir_path(**vars(args))

# Assing arguments
weight = args.weight
pf = args.pf
pa = args.pa
theta_multiplier = args.theta
gamma_multiplier = args.gamma
site_num = args.sitenum
T = args.time
xi = args.xi
dataname = args.dataname
mix_in = args.mix_in
mass_matrix_theta_scale = args.mass_matrix_theta_scale
mass_matrix_gamma_scale = args.mass_matrix_gamma_scale
mass_matrix_weight = args.mass_matrix_weight
symplectic_integrator_num_steps = args.symplectic_integrator_num_steps
stepsize = args.stepsize
scale = args.scale
mode = args.mode

# Solve model with Casadi
casadi_results = solve_with_casadi(
    weight=weight,
    xi=xi,
    pf=pf,
    pa=pa,
    site_num=site_num,
    T=T,
    output_dir=output_dir,
    mix_in=mix_in,
    mass_matrix_theta_scale=mass_matrix_theta_scale,
    mass_matrix_gamma_scale=mass_matrix_gamma_scale,
    mass_matrix_weight=mass_matrix_weight,
    stepsize=stepsize,
    symplectic_integrator_num_steps=symplectic_integrator_num_steps,
    two_param_uncertainty=True,
    scale=scale,
    mode=mode,
)
