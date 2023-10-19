#!/usr/bin/env python
# coding: utf-8

"""
This module provides essential tools to handle data files related to this project.
"""

# Import Required Packages
# ========================

import geopandas as gpd
import pandas as pd
from services.file_service import get_path

# Default Path to the data folder
_DATA_FOLDER = get_path("data", "hmc")


def load_site_data(
    num_sites,
    norm_fac=1e9,
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
    file_path = data_folder / f"hmc_{num_sites}SitesModel.csv"
    df = pd.read_csv(file_path)

    # Extract information
    z_2017 = df[f"z_2017_{num_sites}Sites"].to_numpy()
    zbar_2017 = df[f"zbar_2017_{num_sites}Sites"].to_numpy()
    theta = df[f"theta_{num_sites}Sites"].to_numpy()
    gamma = df[f"gamma_{num_sites}Sites"].to_numpy()
    forestArea_2017_ha = df[f"forestArea_2017_ha_{num_sites}Sites"].to_numpy()

    # Normalize Z data
    zbar_2017 /= norm_fac
    z_2017 /= norm_fac
    forestArea_2017_ha /= norm_fac

    # Read geojson files
    theta_data = gpd.read_file(data_folder / "data_theta.geojson")
    gamma_data = gpd.read_file(data_folder / "data_gamma.geojson")
    id_data = gpd.read_file(data_folder / f"id_{num_sites}.geojson")

    site_theta_2017 = gpd.overlay(id_data, theta_data, how="intersection")
    site_theta_2017_df = site_theta_2017.iloc[:, :-1]
    site_gamma_2017 = gpd.overlay(id_data, gamma_data, how="intersection")
    site_gamma_2017_df = site_gamma_2017.iloc[:, :-1]

    print(f"Data successfully loaded from '{file_path}'")
    return (
        zbar_2017,
        gamma,
        z_2017,
        forestArea_2017_ha,
        theta,
        site_theta_2017_df,
        site_gamma_2017_df,
    )
