import numpy as np
import pandas as pd

from ..services.file_service import get_path


def theta_adj_reg_data(num_sites, theta_df):
    M = 62

    # Get design matrix and its dimensions
    X = theta_df.iloc[:, :8].to_numpy()
    N, K = X.shape

    # Large group indicator
    G = theta_df.iloc[:, 8]
    G = G.astype(int)

    # Municipal level to site level projection matrix
    SG = np.array(
        [(theta_df["id"].to_numpy() == i).astype(int) for i in range(1, num_sites + 1)]
    )

    # Multiply by area overalp weights
    SG = theta_df["muni_site_area"].to_numpy() * SG
    SG = SG / SG.sum(axis=1, keepdims=True)

    return {
        "M_theta": M,
        "X_theta": X,
        "N_theta": N,
        "K_theta": K,
        "SG_theta": SG,
        "G_theta": G,
    }


def gamma_adj_reg_data(num_sites, gamma_df):
    # Get design matrix and its dimensions
    M = 78
    X = gamma_df.iloc[:, :6].to_numpy()
    N, K = X.shape

    # Large group indicator
    G = gamma_df.iloc[:, 7]

    return {
        "M_gamma": M,
        "X_gamma": X,
        "N_gamma": N,
        "K_gamma": K,
        "G_gamma": G,
    }


def baseline_hyperparams(var):
    # Load baseline hyperparameters
    data_dir = get_path("data", "calibration", "hmc")
    df = pd.read_csv(data_dir / f"{var}_posterior_means_and_variances.csv")

    beta_mean = df[df["Parameter"].str.startswith("beta_")]["Posterior Mean"].values
    V_mean = df[df["Parameter"].str.startswith("V_")]["Posterior Mean"].values

    # Construct coefficients vcov
    beta_cov = np.zeros((len(beta_mean), len(beta_mean)))
    for i in range(len(beta_mean)):
        for j in range(len(beta_mean)):
            param_name = f"var_beta_{i+1}_{j+1}"
            beta_cov[i, j] = df[df["Parameter"] == param_name][
                "Posterior Variance"
            ].values[0]

    V_var = df[df["Parameter"].str.startswith("V_")]["Posterior Variance"].values

    return {
        f"beta_{var}_mean": beta_mean,
        f"beta_{var}_vcov": beta_cov,
        f"V_{var}_mean": V_mean,
        f"V_{var}_var": V_var,
    }




