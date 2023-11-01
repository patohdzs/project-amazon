import numpy as np


def theta_adj_reg_data(num_sites, theta_df):
    # Filter out null values
    theta_df = theta_df[theta_df["zbar_2017_muni"].notna()]

    # Get design matrix and its dimensions
    X = theta_df.iloc[:, 1:9].to_numpy()
    N, K = X.shape

    # Get weighted grouped average matrix
    G = np.array(
        [(theta_df["id"].to_numpy() == i).astype(int) for i in range(1, num_sites + 1)]
    )
    G = theta_df["zbar_2017_muni"].to_numpy() * G
    G = G / G.sum(axis=1, keepdims=True)

    return {
        "X_theta": X,
        "N_theta": N,
        "K_theta": K,
        "G_theta": G,
    }


def gamma_adj_reg_data(num_sites, gamma_df):
    # Get design matrix and its dimensions
    X = gamma_df.iloc[:, 1:6].to_numpy()
    N, K = X.shape

    # Get grouped average matrix
    G = np.array(
        [(gamma_df["id"].to_numpy() == i).astype(int) for i in range(1, num_sites + 1)]
    )
    G = G / G.sum(axis=1, keepdims=True)

    return {
        "X_gamma": X,
        "N_gamma": N,
        "K_gamma": K,
        "G_gamma": G,
    }
