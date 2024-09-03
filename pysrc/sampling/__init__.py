import numpy as np
import pandas as pd
from scipy.linalg import cholesky
from ..services.file_service import get_path

def theta_adj_reg_data(num_sites, theta_df):
    # Get design matrix and its dimensions
    X = theta_df.iloc[:, :8].to_numpy()
    N, K = X.shape

    # Get site indicator matrix
    G = np.array(
        [(theta_df["id"].to_numpy() == i).astype(int) for i in range(1, num_sites + 1)]
    )

    # Multiply by area overalp weights
    G = theta_df["muni_site_area"].to_numpy() * G
    G = G / G.sum(axis=1, keepdims=True)
    return {
        "X_theta": X,
        "N_theta": N,
        "K_theta": K,
        "G_theta": G,
    }


def gamma_adj_reg_data(num_sites, gamma_df):
    # Get design matrix and its dimensions
    X = gamma_df.iloc[:, :5].to_numpy()
    N, K = X.shape

    # Get grouped average matrix
    G = np.array(
        [(gamma_df["id"].to_numpy() == i).astype(int) for i in range(1, num_sites + 1)]
    )

    # Multiply by area overalp weights
    G = gamma_df["muni_site_area"].to_numpy() * G
    G = G / G.sum(axis=1, keepdims=True)

    return {
        "X_gamma": X,
        "N_gamma": N,
        "K_gamma": K,
        "G_gamma": G,
    }



def gibbs_sampling():
    # Load the parameters from the CSV file
    data_dir = get_path("data", "calibration", "hmc")
    df = pd.read_csv(data_dir / 'posterior_means_and_variances.csv')

    # Extract the posterior means
    beta_mean = df[df['Parameter'].str.startswith('beta_') & df['Posterior Mean'].notna()]['Posterior Mean'].values
    eta_mean = df[df['Parameter'] == 'eta']['Posterior Mean'].values[0]
    nu_mean = df[df['Parameter'] == 'nu']['Posterior Mean'].values[0]
    V_mean = df[df['Parameter'].str.startswith('V_') & df['Posterior Mean'].notna()]['Posterior Mean'].values

    # Reconstruct the beta variance (covariance matrix)
    beta_cov = np.zeros((len(beta_mean), len(beta_mean)))
    for i in range(len(beta_mean)):
        for j in range(len(beta_mean)):
            param_name = f'var_beta_{i+1}_{j+1}'
            beta_cov[i, j] = df[df['Parameter'] == param_name]['Posterior Variance'].values[0]

    # Extract the posterior variances
    eta_var = df[df['Parameter'] == 'eta']['Posterior Variance'].values[0]
    nu_var = df[df['Parameter'] == 'nu']['Posterior Variance'].values[0]
    V_var = df[df['Parameter'].str.startswith('V_') & df['Posterior Variance'].notna()]['Posterior Variance'].values


    return_dict = {
    "a_gamma": (eta_mean ** 2) / eta_var,
    "b_gamma": eta_var / eta_mean,
    "a_nu":(nu_mean ** 2) / nu_var,
    "b_nu":nu_var / nu_mean,
    "gamma_mean":beta_mean,
    "gamma_vcov":beta_cov,
    "V_mean":V_mean,
    "V_var":V_var,
    }
    
    return return_dict