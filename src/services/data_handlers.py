#!/usr/bin/env python
# coding: utf-8

"""
This module provides essential tools to handle data files related to this project.
"""

# Import Required Packages
# ========================
import os
import pandas as pd


# Default Path to the data folder
_DATA_FOLDER = os.path.join(os.path.dirname(__file__), "calibration")


def load_site_data(
    site_num,
    norm_fac=1.0,
    data_folder=_DATA_FOLDER,
):
    """
    Load site data given the number of sites `site_num`.
    The normalization factor is used to scale `z_2017` and `zbar_2017` (multiplicative) data.

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
    file = os.path.join(data_folder, f"calibration_{n}SitesModel.csv")
    df = pd.read_csv(file)

    # Extract information
    z_2017 = df[f"z_2017_{n}Sites"].to_numpy()
    zbar_2017 = df[f"zbar_2017_{n}Sites"].to_numpy()
    gamma = df[f"gamma_{n}Sites"].to_numpy()
    gammaSD = df[f"gammaSD_{n}Sites"].to_numpy()
    forestArea_2017_ha = df[f"forestArea_2017_ha_{n}Sites"].to_numpy()
    theta = df[f"theta_{n}Sites"].to_numpy()
    try:
        # Old files do not provide thetaSD
        thetaSD = df[f"thetaSD_{n}Sites"].to_numpy()
    except KeyError:
        thetaSD = None

    # Normalize Z data
    zbar_2017 /= norm_fac
    z_2017 /= norm_fac

    print(f"Data successfully loaded from '{file}'")
    return (
        zbar_2017,
        gamma,
        gammaSD,
        z_2017,
        forestArea_2017_ha,
        theta,
        thetaSD,
    )
