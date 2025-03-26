import pandas as pd
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
from scipy.stats import gaussian_kde
from scipy.stats import entropy
import geopandas as gpd
import pandas as pd

from matplotlib.backends.backend_pdf import PdfPages
import sys
from pysrc.services.file_service import get_path

import argparse
parser = argparse.ArgumentParser(description="shadow price calculation")
parser.add_argument("--pee",type=float,default=5)
parser.add_argument("--xi",type=float,default=5)
parser.add_argument("--sites",type=int,default=78)

args = parser.parse_args()
pee=args.pee
xi=args.xi
num_sites=args.sites

solver="gurobi"
pa=41.11



result_folder = os.path.join(
    str(get_path("output")),
    "sampling",
    solver,
    f"{num_sites}sites",
    f"pa_{pa}",
    f"xi_{xi}",
)
prior_folder = os.path.join(
    str(get_path("output")),
    "sampling",
    solver,
    f"{num_sites}sites",
    f"pa_{pa}",
    "xi_10000",
)

with open(result_folder + f"/pe_{pee}/results.pcl", "rb") as f:
    b0 = pickle.load(f)

with open(result_folder + f"/pe_{pee+15}/results.pcl", "rb") as f:
    b15 = pickle.load(f)

with open(prior_folder + f"/pe_{pee+15}/results.pcl", "rb") as f:
    results_unadjusted = pickle.load(f)
    
    
    
theta_unadjusted = results_unadjusted["final_sample"][:, :num_sites]
gamma_unadjusted = results_unadjusted["final_sample"][:, num_sites:]
theta_adjusted_b0 = b0["final_sample"][:, :num_sites]
gamma_adjusted_b0 = b0["final_sample"][:, num_sites:]
theta_adjusted_b15 = b15["final_sample"][:, :num_sites]
gamma_adjusted_b15 = b15["final_sample"][:, num_sites:]




b = [0, 15]

kl_divergences = []


theta_list = [theta_adjusted_b0, theta_adjusted_b15]  
gamma_list = [gamma_adjusted_b0, gamma_adjusted_b15]  


for order, (theta_hmc, gamma_hmc) in enumerate(zip(theta_list, gamma_list)):
    for idx in range(theta_hmc.shape[1]):  
        
        
        joint_samples_unadjusted = np.vstack([theta_unadjusted[:, idx], gamma_unadjusted[:, idx]]).T  # shape (16000, 2)
        joint_samples_adjusted = np.vstack([theta_hmc[:, idx], gamma_hmc[:, idx]]).T  # shape (16000, 2)

        kde_func_unadjusted = gaussian_kde(joint_samples_unadjusted.T, bw_method='scott')
        kde_func_adjusted = gaussian_kde(joint_samples_adjusted.T, bw_method='scott')

        theta_min = min(theta_unadjusted[:, idx].min(), theta_hmc[:, idx].min())
        theta_max = max(theta_unadjusted[:, idx].max(), theta_hmc[:, idx].max())
        gamma_min = min(gamma_unadjusted[:, idx].min(), gamma_hmc[:, idx].min())
        gamma_max = max(gamma_unadjusted[:, idx].max(), gamma_hmc[:, idx].max())

        # Define a grid of points in the 2D (theta, gamma) space
        theta_grid = np.linspace(theta_min, theta_max, 100)
        gamma_grid = np.linspace(gamma_min, gamma_max, 100)
        theta_mesh, gamma_mesh = np.meshgrid(theta_grid, gamma_grid)
        grid_points = np.vstack([theta_mesh.ravel(), gamma_mesh.ravel()])  # shape (2, 10000)

        density_unadjusted = kde_func_unadjusted(grid_points)
        density_adjusted = kde_func_adjusted(grid_points)

        density_unadjusted += 1e-20
        density_adjusted += 1e-20


        kl_div = entropy(density_unadjusted, density_adjusted)
        print(f'Site {idx + 1} Joint (Theta, Gamma) KL Divergence (b = {b[order]}): {kl_div}')

        kl_divergences.append({'Parameter': f'Joint_theta_gamma_b{b[order]}', 'Site': idx + 1, 'KL_Divergence': kl_div})



output_folder = str(get_path("output")) + f"/figures/density/site_{num_sites}/xi{xi}/"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)


kl_df = pd.DataFrame(kl_divergences)
kl_df.to_csv(output_folder+'kl_divergences_theta_gamma.csv', index=False)


top_kl_divergences = []

for parameter in kl_df['Parameter'].unique():
    top_sites = kl_df[kl_df['Parameter'] == parameter].nlargest(2, 'KL_Divergence')
    top_kl_divergences.append(top_sites)

# Concatenate all the top results
print("top2 kl",top_kl_divergences)
