#!/usr/bin/env python

# Import Required Packages
# ========================
import argparse
import os
import pickle
import sys
from pathlib import Path

import plots
import seaborn as sns
from services.data_service import load_site_data

sys.path.append(os.path.abspath("src"))
sns.set(font_scale=1.2)

parser = argparse.ArgumentParser(description="parameter settings")
parser.add_argument("--weight", type=float, default=0.25)
parser.add_argument("--xi", type=float, default=0.01)
parser.add_argument("--pf", type=float, default=20.76)
parser.add_argument("--pa", type=float, default=44.75)
parser.add_argument("--theta", type=float, default=1.0)
parser.add_argument("--gamma", type=float, default=1.0)
parser.add_argument("--sitenum", type=int, default=10)
parser.add_argument("--time", type=int, default=200)
parser.add_argument("--dataname", type=str, default="tests")
parser.add_argument("--mix_in", type=int, default=2)
parser.add_argument("--mass_matrix_theta_scale", type=float, default=1.0)
parser.add_argument("--mass_matrix_gamma_scale", type=float, default=1.0)
parser.add_argument("--mass_matrix_weight", type=float, default=0.1)
parser.add_argument("--symplectic_integrator_num_steps", type=int, default=2)
parser.add_argument("--stepsize", type=float, default=0.1)
parser.add_argument("--scale", type=float, default=1.0)
parser.add_argument("--mode", type=float, default=1.0)


args = parser.parse_args()
weight = args.weight
pf = args.pf
pa = args.pa
theta_multiplier = args.theta
gamma_multiplier = args.gamma
sitenum = args.sitenum
time = args.time
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

workdir = os.getcwd()
output_dir = (
    workdir
    + "/output/"
    + dataname
    + "/scale_"
    + str(scale)
    + "_mode_"
    + str(mode)
    + "/pf_"
    + str(pf)
    + "_pa_"
    + str(pa)
    + "_time_"
    + str(time)
    + "/theta_"
    + str(theta_multiplier)
    + "_gamma_"
    + str(gamma_multiplier)
    + "/sitenum_"
    + str(sitenum)
    + "_xi_"
    + str(xi)
    + "/mix_in_"
    + str(mix_in)
    + "_mm_theta_scale_"
    + str(mass_matrix_theta_scale)
    + "_mm_gamma_scale_"
    + str(mass_matrix_gamma_scale)
    + "_num_steps_"
    + str(symplectic_integrator_num_steps)
    + "_stepsize_"
    + str(stepsize)
    + "/weight_"
    + str(weight)
    + "_mass_matrix_weight_"
    + str(mass_matrix_weight)
    + "/"
)
plotdir = (
    workdir
    + "/plot/"
    + dataname
    + "/scale_"
    + str(scale)
    + "_mode_"
    + str(mode)
    + "/pf_"
    + str(pf)
    + "_pa_"
    + str(pa)
    + "_time_"
    + str(time)
    + "/theta_"
    + str(theta_multiplier)
    + "_gamma_"
    + str(gamma_multiplier)
    + "/sitenum_"
    + str(sitenum)
    + "_xi_"
    + str(xi)
    + "/mix_in_"
    + str(mix_in)
    + "_mm_theta_scale_"
    + str(mass_matrix_theta_scale)
    + "_mm_gamma_scale_"
    + str(mass_matrix_gamma_scale)
    + "_num_steps_"
    + str(symplectic_integrator_num_steps)
    + "_stepsize_"
    + str(stepsize)
    + "/weight_"
    + str(weight)
    + "_mass_matrix_weight_"
    + str(mass_matrix_weight)
    + "/"
)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
if not os.path.exists(plotdir):
    os.makedirs(plotdir)


(
    zbar_2017,
    gamma,
    z_2017,
    forestArea_2017_ha,
    theta,
    gamma_coe,
    gamma_coe_sd,
    theta_coe,
    theta_coe_sd,
    gamma_vcov_array,
    theta_vcov_array,
    site_theta_2017_df,
    site_gamma_2017_df,
) = load_site_data(
    sitenum,
    norm_fac=1e9,
)


with open(output_dir + "results.pcl", "rb") as f:
    # Load the data from the file
    results = pickle.load(f)


# Plot absolute and percentage error
plots.traceplot_abs_error(results, Path(plotdir))
plots.traceplot_pct_error(results, Path(plotdir))


# Plot theta and gamma convergence
plots.traceplot_params_pct_error(results, Path(plotdir))
plots.traceplot_gamma(results, Path(plotdir))
plots.traceplot_theta(results, Path(plotdir))
