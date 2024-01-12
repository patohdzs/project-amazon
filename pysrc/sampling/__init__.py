import numpy as np


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

    # Get inv vector of weights
    w = 1 / theta_df["weights"]

    return {
        "X_theta": X,
        "N_theta": N,
        "K_theta": K,
        "G_theta": G,
        "w_theta": w,
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
