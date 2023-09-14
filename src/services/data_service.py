#!/usr/bin/env python
# coding: utf-8

"""
This module provides essential tools to handle data files related to this project.
"""

# Import Required Packages
# ========================
import os
import pandas as pd
from .file_service import get_path


def load_site_data(site_num, norm_fac=1.0):
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
        - thetaSD:
    """
    # Read data file
    n = site_num
    file_path = get_path(
        "data", "calibration", f"{n}SitesModel", f"calibration_{n}SitesModel.csv"
    )
    df = pd.read_csv(file_path)

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

    print(f"Data successfully loaded from '{file_path}'")
    return (
        zbar_2017,
        gamma,
        gammaSD,
        z_2017,
        forestArea_2017_ha,
        theta,
        thetaSD,
    )
