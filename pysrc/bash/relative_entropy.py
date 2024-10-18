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

solver="gams"
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
    
    
    
theta_unadjusted = results_unadjusted["final_sample"][:16000, :num_sites]
gamma_unadjusted = results_unadjusted["final_sample"][:16000, num_sites:]
theta_adjusted_b0 = b0["final_sample"][:16000, :num_sites]
gamma_adjusted_b0 = b0["final_sample"][:16000, num_sites:]
theta_adjusted_b15 = b15["final_sample"][:16000, :num_sites]
gamma_adjusted_b15 = b15["final_sample"][:16000, num_sites:]




b = [0, 15]
theta_list = [theta_adjusted_b0, theta_adjusted_b15]
gamma_list = [gamma_adjusted_b0, gamma_adjusted_b15]
theta_hmc_list = ['theta_b0', 'theta_b15']
gamma_hmc_list = ['gamma_b0', 'gamma_b15']


kl_divergences = []

for order, theta_hmc in enumerate(theta_list):
    for idx in range(theta_hmc.shape[1]):
        kde_func_unadjusted = gaussian_kde(theta_unadjusted[:, idx], bw_method='scott')
        kde_func_adjusted = gaussian_kde(theta_hmc[:, idx], bw_method='scott')

        common_grid = np.linspace(min(theta_unadjusted[:, idx].min(), theta_hmc[:, idx].min()),
                                  max(theta_unadjusted[:, idx].max(), theta_hmc[:, idx].max()),
                                  1000)
        density_unadjusted = kde_func_unadjusted(common_grid)
        density_adjusted = kde_func_adjusted(common_grid)
        density_unadjusted += 1e-20
        density_adjusted += 1e-20

        kl_div = entropy(density_unadjusted, density_adjusted)
        print(f'Site {idx + 1} Theta KL Divergence (b = {b[order]}): {kl_div}')

        kl_divergences.append({'Parameter': theta_hmc_list[order], 'Site': idx + 1, 'KL_Divergence': kl_div})


for order, gamma_hmc in enumerate(gamma_list):
    for idx in range(gamma_hmc.shape[1]):
        kde_func_unadjusted = gaussian_kde(gamma_unadjusted[:, idx], bw_method='scott')
        kde_func_adjusted = gaussian_kde(gamma_hmc[:, idx], bw_method='scott')

        common_grid = np.linspace(min(gamma_unadjusted[:, idx].min(), gamma_hmc[:, idx].min()),
                                  max(gamma_unadjusted[:, idx].max(), gamma_hmc[:, idx].max()),
                                  1000)
        density_unadjusted = kde_func_unadjusted(common_grid)
        density_adjusted = kde_func_adjusted(common_grid)
        density_unadjusted += 1e-20
        density_adjusted += 1e-20

        kl_div = entropy(density_unadjusted, density_adjusted)
        print(f'Site {idx + 1} Gamma KL Divergence (b = {b[order]}): {kl_div}')

        kl_divergences.append({'Parameter': gamma_hmc_list[order], 'Site': idx + 1, 'KL_Divergence': kl_div})



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
