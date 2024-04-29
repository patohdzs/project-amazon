import argparse
import os

from pysrc.sampling import mpc

# This script will sample price paths for constrained and unconstrained models

# Read arguments from stdin
parser = argparse.ArgumentParser(description="parameter settings")
parser.add_argument("--sitenum", type=int, default=78)
# Parse arguments
args = parser.parse_args()


# Define the root folders
root_folder = os.getcwd() + "/gams_file/mpc"

computation_folder = os.path.join(root_folder, "computation")
simulation_folder = os.path.join(root_folder, "simulation")
calculation_folder_unconstrained = os.path.join(
    root_folder + "/computation", "model_unconstrained"
)
mc_samples_folder_unconstrained = os.path.join(
    root_folder + "/simulation", "sample_unconstrained"
)
calculation_folder_constrained = os.path.join(
    root_folder + "/computation", "model_constrained"
)
mc_samples_folder_constrained = os.path.join(
    root_folder + "/simulation", "sample_constrained"
)


os.makedirs(computation_folder, exist_ok=True)
os.makedirs(simulation_folder, exist_ok=True)
os.makedirs(calculation_folder_unconstrained, exist_ok=True)
os.makedirs(mc_samples_folder_unconstrained, exist_ok=True)
os.makedirs(calculation_folder_constrained, exist_ok=True)
os.makedirs(mc_samples_folder_constrained, exist_ok=True)

# get mc samples
mpc.sample_price_paths(location=mc_samples_folder_unconstrained, var="uncon")
mpc.sample_price_paths(location=mc_samples_folder_constrained, var="con")
print("sampling is done")


# get gdx files
mpc.gdx_files(num_sites=args.sitenum, location=calculation_folder_unconstrained)
mpc.gdx_files(num_sites=args.sitenum, location=calculation_folder_constrained)

print("gdx file is done")


# move samples into computation folders
mpc.paste_file(
    root=root_folder,
    ori=mc_samples_folder_unconstrained,
    des=calculation_folder_unconstrained,
)
mpc.paste_file(
    root=root_folder,
    ori=mc_samples_folder_constrained,
    des=calculation_folder_constrained,
    model="constrained",
)


print("file movement is done")
