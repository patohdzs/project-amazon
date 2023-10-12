#!/usr/bin/env python
# coding: utf-8

"""
This module provides essential tools to handle data files related to this project.
"""

# Import Required Packages
# ========================

import geopandas as gpd
import numpy as np
import pandas as pd

from .file_service import get_path

# Default Path to the data folder
_DATA_FOLDER = get_path("data", "hmc")


def load_site_data(
    site_num,
    norm_fac=1.0,
    data_folder=_DATA_FOLDER,
):
    """
    Load site data given the number of sites `site_num`.
    The normalization factor is used to scale `z_2017`
    and `zbar_2017` (multiplicative) data.

    :returns:
        - zbar_2017:
        - gamma:
        - gammaSD:
        - z_2017:
        - forestArea_2017_ha:
        - theta:
    """
    # Read data file
    n = site_num
    file_path = data_folder / f"hmc_{n}SitesModel.csv"
    df = pd.read_csv(file_path)

    # Extract information
    z_2017 = df[f"z_2017_{n}Sites"].to_numpy()
    zbar_2017 = df[f"zbar_2017_{n}Sites"].to_numpy()
    gamma = df[f"gamma_{n}Sites"].to_numpy()
    # gammaSD            = df[f'gammaSD_{n}Sites'].to_numpy()
    forestArea_2017_ha = df[f"forestArea_2017_ha_{n}Sites"].to_numpy()
    theta = df[f"theta_{n}Sites"].to_numpy()

    # Normalize Z data
    zbar_2017 /= norm_fac
    z_2017 /= norm_fac

    # Read in prior samples
    gamma_coe_prior = pd.read_csv(data_folder / "gamma_coe.csv").to_numpy()
    theta_coe_prior = pd.read_csv(data_folder / "theta_coe.csv").to_numpy()

    # Compute means
    gamma_coe_mean = gamma_coe_prior.mean(axis=0)
    theta_coe_mean = theta_coe_prior.mean(axis=0)

    # Compute covariances
    gamma_coe_vcov = np.cov(gamma_coe_prior, rowvar=False)
    theta_coe_vcov = np.cov(theta_coe_prior, rowvar=False)

    # Read geojson files
    data_theta = gpd.read_file(data_folder / "data_theta.geojson")
    data_gamma = gpd.read_file(data_folder / "data_gamma.geojson")

    file_path_id = data_folder / f"id_{n}.geojson"
    data_id = gpd.read_file(file_path_id)

    site_theta_2017 = gpd.overlay(data_id, data_theta, how="intersection")
    site_theta_2017_df = site_theta_2017.iloc[:, :-1]
    site_gamma_2017 = gpd.overlay(data_id, data_gamma, how="intersection")
    site_gamma_2017_df = site_gamma_2017.iloc[:, :-1]

    print(f"Data successfully loaded from '{file_path}'")
    return (
        zbar_2017,
        gamma,
        z_2017,
        forestArea_2017_ha,
        theta,
        gamma_coe_mean,
        theta_coe_mean,
        gamma_coe_vcov,
        theta_coe_vcov,
        site_theta_2017_df,
        site_gamma_2017_df,
    )


def load_coef_prior_samples():
    # Load coef prior samples
    beta_theta_prior_samples = pd.read_csv(
        get_path("data", "hmc", "theta_coe.csv")
    ).to_numpy()

    beta_gamma_prior_samples = pd.read_csv(
        get_path("data", "hmc", "gamma_coe.csv")
    ).to_numpy()

    return beta_theta_prior_samples, beta_gamma_prior_samples
