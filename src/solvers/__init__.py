#!/usr/bin/env python
# coding: utf-8

"""
This module presents solution approaches for solving the bilevel
optimization problem with MCMC sampling for the inner optimization problem.

The outer optimization problem is either solved with CASADI or GAMS.
"""


# Import Required Packages
# ========================
import os
import sys

import numpy as np

# MCMC (HMC) sampling routines
mcmc_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "src/mcmc")
if mcmc_path not in sys.path:
    sys.path.append(mcmc_path)


# Data Hanlder (.data_handlers.load_site_data)
sys.path.append(os.path.abspath("src"))


# Usefule Variables:
_DEBUG = False


def theta_fitted(theta_coe, theta_dataframe):
    # Copy df
    theta_data = theta_dataframe.copy()

    # Compute fitted values
    theta_data["fitted_value"] = np.exp(
        (theta_data.iloc[:, 1:9] * theta_coe).sum(axis=1)
    )

    # Subset and filter
    theta_data = theta_data[["id", "zbar_2017_muni", "fitted_value"]]
    theta_data = theta_data[theta_data["zbar_2017_muni"].notna()]

    # Weighted average function
    aux_price_2017 = 44.9736197781184

    def weighted_mean(group):
        return np.average(
            group["fitted_value"] / aux_price_2017, weights=group["zbar_2017_muni"]
        )

    # Take weighted average
    result = (
        theta_data.groupby("id")
        .apply(weighted_mean)
        .reset_index(name="theta2017_Sites")
    )

    # Return as numpy array
    return result["theta2017_Sites"].to_numpy()


def gamma_fitted(gamma_coe, gamma_dataframe):
    # Copy df
    gamma_data = gamma_dataframe.copy()

    # Compute fitted values
    gamma_data["fitted_value"] = np.exp(
        (gamma_data.iloc[:, 1:6] * gamma_coe).sum(axis=1)
    )

    # Subset columns
    gamma_data = gamma_data[["id", "fitted_value"]]

    # Group by id and compute weighted mean for each group
    result = (
        gamma_data.groupby("id")["fitted_value"]
        .mean()
        .reset_index(name="gamma2017_Sites")
    )
    return result["gamma2017_Sites"].to_numpy()


def log_density_function(
    uncertain_vals,
    uncertain_vals_mean,
    block_matrix,
    N,  # Number of control intervals
    alpha,  # Mean reversion coefficient
    sol_val_X,
    sol_val_Ua,  # Squared total control adjustments; dimensions T x 1
    sol_val_Up,  # U control; dimensions T x I
    zbar_2017,
    forestArea_2017_ha,  # Forest area in 2017
    norm_fac,  # Normalization factor
    alpha_p_Adym,  # Vector of alpha^Adym
    Bdym,
    leng,
    T,  # Time horizon
    ds_vect,  # Time discounting vector
    zeta,  # Adjustment costs parameter
    xi,  # Penalty on KL divergence
    kappa,  # Effect of cattle farming on emissions
    pa,  # Price of cattle output
    pf,  # Shadow price of carbon emissions
    site_theta_2017_df,
    site_gamma_2017_df,
):
    """
    Define a function to evaluate log-density of the objective/posterior distribution.

    Some of the input parameters are updated at each cycle of the outer loop
    (optimization loop),and it becomes then easier/cheaper to udpate the function stamp
    and keep it separate here

    Note that the log-density is the logarithm of the target density discarding any
    normalization factor
    """
    # Flatten time-discounting vector
    ds_vect = np.asarray(ds_vect).flatten()

    # Flatten uncertain values
    uncertain_vals = np.asarray(uncertain_vals).flatten()

    # Unpacking uncertain values
    theta_coe_vals = uncertain_vals[:8]
    gamma_coe_vals = uncertain_vals[8:]

    # Computing fitted values for theta and gamma
    theta_fit = theta_fitted(
        theta_coe=theta_coe_vals, theta_dataframe=site_theta_2017_df
    )
    gamma_fit = gamma_fitted(
        gamma_coe=gamma_coe_vals, gamma_dataframe=site_gamma_2017_df
    )

    # Num of theta_fitted values
    size = theta_fit.size

    # Full fitted values vector (theta, gamma)
    beta_vals = np.concatenate((theta_fit, gamma_fit))

    # Check for complex-valued uncertain_vals
    if np.iscomplexobj(beta_vals):
        raise TypeError(f"beta_vals contains complex numbers {beta_vals}")

    # Carbon captured at time zero
    x0_vals = beta_vals[size:].T.dot(forestArea_2017_ha) / norm_fac
    X_zero = np.sum(x0_vals) * np.ones(leng)

    shifted_X = sol_val_X[0:size, :-1].copy()
    for j in range(N):
        shifted_X[:, j] = zbar_2017 - shifted_X[:, j]
    omega = np.dot(beta_vals[size:], alpha * shifted_X - sol_val_Up)

    X_dym = np.zeros(T + 1)
    X_dym[0] = np.sum(x0_vals)
    X_dym[1:] = alpha_p_Adym * X_zero + np.dot(Bdym, omega.T)

    z_shifted_X = sol_val_X[0:size, :].copy()
    scl = pa * beta_vals[:size] - pf * kappa

    for j in range(N + 1):
        z_shifted_X[:, j] *= scl

    # Adjustment costs
    term_1 = -np.sum(ds_vect[0:T] * sol_val_Ua) * zeta / 2

    # Value of emissions absorbed
    term_2 = np.sum(ds_vect[0:T] * (X_dym[1:] - X_dym[0:-1])) * pf

    # Value of cattle output minus cost of emissions
    term_3 = np.sum(ds_vect * np.sum(z_shifted_X, axis=0))

    # Overall objective value
    obj_val = term_1 + term_2 + term_3

    # Computing log prior density term
    log_prior_density = _log_normal_prior(
        uncertain_vals, uncertain_vals_mean, block_matrix
    )

    # Computing potential energy
    log_density_val = -1.0 / xi * obj_val + log_prior_density
    log_density_val = float(log_density_val)

    return log_density_val


def _log_normal_prior(x, mu, Sigma):
    # Log prior density assuming x (prior)~ N(mu, Sigma)
    mean_zero_x = x - mu
    return -0.5 * np.dot(mean_zero_x, np.linalg.inv(Sigma).dot(mean_zero_x))
