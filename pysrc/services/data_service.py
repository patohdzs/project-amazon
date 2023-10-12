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

    # mean and sd for coefficients

    gamma_coe = np.array(
        [
            df["gamma_cons"].to_numpy()[0],
            df["gamma_log_precip"].to_numpy()[0],
            df["gamma_log_temp"].to_numpy()[0],
            df["gamma_loglat"].to_numpy()[0],
            df["gamma_loglon"].to_numpy()[0],
        ]
    )
    gamma_coe_sd = np.array(
        [
            df["gamma_sd_cons"].to_numpy()[0],
            df["gamma_sd_log_precip"].to_numpy()[0],
            df["gamma_sd_log_temp"].to_numpy()[0],
            df["gamma_sd_loglat"].to_numpy()[0],
            df["gamma_sd_loglon"].to_numpy()[0],
        ]
    )
    theta_coe = np.array(
        [
            df["theta_cons"].to_numpy()[0],
            df["theta_precip"].to_numpy()[0],
            df["theta_temp"].to_numpy()[0],
            df["theta_temp2"].to_numpy()[0],
            df["theta_lat"].to_numpy()[0],
            df["theta_lat2"].to_numpy()[0],
            df["theta_gateprice"].to_numpy()[0],
            df["theta_distance"].to_numpy()[0],
        ]
    )
    theta_coe_sd = np.array(
        [
            df["theta_sd_cons"].to_numpy()[0],
            df["theta_sd_precip"].to_numpy()[0],
            df["theta_sd_temp"].to_numpy()[0],
            df["theta_sd_temp2"].to_numpy()[0],
            df["theta_sd_lat"].to_numpy()[0],
            df["theta_sd_lat2"].to_numpy()[0],
            df["theta_sd_gateprice"].to_numpy()[0],
            df["theta_sd_distance"].to_numpy()[0],
        ]
    )

    # Normalize Z data
    zbar_2017 /= norm_fac
    z_2017 /= norm_fac

    file_path_2 = data_folder / "gamma_vcov.csv"
    df2 = pd.read_csv(file_path_2, header=None)

    file_path_3 = data_folder / "theta_vcov.csv"
    df3 = pd.read_csv(file_path_3, header=None)

    gamma_vcov_array = df2.values
    theta_vcov_array = df3.values

    file_path_data_theta = data_folder / "data_theta.geojson"
    data_theta = gpd.read_file(file_path_data_theta)

    file_path_data_gamma = data_folder / "data_gamma.geojson"
    data_gamma = gpd.read_file(file_path_data_gamma)

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
        gamma_coe,
        gamma_coe_sd,
        theta_coe,
        theta_coe_sd,
        gamma_vcov_array,
        theta_vcov_array,
        site_theta_2017_df,
        site_gamma_2017_df,
    )