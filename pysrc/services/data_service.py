import geopandas as gpd
import pandas as pd

from ..services.file_service import get_path

# Default Path to the data folder
_DATA_FOLDER = get_path("data", "hmc")


def load_site_data(
    num_sites,
    norm_fac=1e9,
    data_folder=_DATA_FOLDER,
):
    # Read data file
    file_path = data_folder / f"hmc_{num_sites}SitesModel.csv"
    df = pd.read_csv(file_path)

    # Extract information
    z_2017 = df[f"z_2017_{num_sites}Sites"].to_numpy()
    zbar_2017 = df[f"zbar_2017_{num_sites}Sites"].to_numpy()
    theta = df[f"theta_{num_sites}Sites"].to_numpy()
    gamma = df[f"gamma_{num_sites}Sites"].to_numpy()
    forest_area_2017 = df[f"forestArea_2017_ha_{num_sites}Sites"].to_numpy()

    # Normalize Z data
    zbar_2017 /= norm_fac
    z_2017 /= norm_fac
    forest_area_2017 /= norm_fac

    # Read municipal level data
    municipal_theta_df = gpd.read_file(data_folder / "muni_data_theta.geojson")
    municipal_gamma_df = gpd.read_file(data_folder / "muni_data_gamma.geojson")

    # Read site level data
    site_theta_2017 = gpd.read_file(
        data_folder / f"site_{num_sites}_data_theta.geojson"
    )
    site_gamma_2017 = gpd.read_file(
        data_folder / f"site_{num_sites}_data_gamma.geojson"
    )

    # Remove geometries
    site_theta_2017_df = site_theta_2017.iloc[:, :-1]
    site_gamma_2017_df = site_gamma_2017.iloc[:, :-1]

    print(f"Data successfully loaded from {data_folder}")
    return (
        zbar_2017,
        gamma,
        z_2017,
        forest_area_2017,
        theta,
        site_theta_2017_df,
        site_gamma_2017_df,
        municipal_theta_df,
        municipal_gamma_df,
    )
