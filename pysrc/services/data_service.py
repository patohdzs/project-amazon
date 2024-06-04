import geopandas as gpd
import numpy as np
import pandas as pd



from ..services.file_service import get_path


def load_site_data(num_sites: int, norm_fac: float = 1e9):
    # Set data directory
    data_dir = get_path("data", "hmc")

    # Read data file
    file_path = data_dir / f"hmc_{num_sites}SitesModel.csv"
    df = pd.read_csv(file_path)

    # Extract information
    z_2017 = df[f"z_2017_{num_sites}Sites"].to_numpy()
    zbar_2017 = df[f"zbar_2017_{num_sites}Sites"].to_numpy()
    forest_area_2017 = df[f"forestArea_2017_ha_{num_sites}Sites"].to_numpy()

    # Normalize Z and forest data
    zbar_2017 /= norm_fac
    z_2017 /= norm_fac
    forest_area_2017 /= norm_fac

    # Read municipal level data
    municipal_theta_df = gpd.read_file(data_dir / "muni_data_theta.geojson")
    municipal_gamma_df = gpd.read_file(data_dir / "muni_data_gamma.geojson")

    # Read site level data
    site_theta_2017 = gpd.read_file(data_dir / f"site_{num_sites}_data_theta.geojson")
    site_gamma_2017 = gpd.read_file(data_dir / f"site_{num_sites}_data_gamma.geojson")

    # Remove geometries
    site_theta_2017_df = site_theta_2017.iloc[:, :-1]
    site_gamma_2017_df = site_gamma_2017.iloc[:, :-1]

    print(f"Data successfully loaded from {data_dir}")
    return (
        zbar_2017,
        z_2017,
        forest_area_2017,
        site_theta_2017_df,
        site_gamma_2017_df,
        municipal_theta_df,
        municipal_gamma_df,
    )


def load_price_data():
    # Read data file
    file_path = get_path("data", "hmc") / "seriesPriceCattle_prepared.csv"
    df = pd.read_csv(file_path)
    average_prices = df.groupby("year")["price_real_mon_cattle"].mean()
    p_a_list = np.array(average_prices)
    return p_a_list


def load_site_data_1995(num_sites: int, norm_fac: float = 1e9):
    from pysrc.sampling import baseline
    # Set data directory
    data_dir = get_path("data", "hmc")

    # Read data file
    file_path = data_dir / f"hmc_{num_sites}SitesModel.csv"
    df = pd.read_csv(file_path)

    # Extract information
    z_1995 = df[f"z_1995_{num_sites}Sites"].to_numpy()
    z_2008 = df[f"z_2008_{num_sites}Sites"].to_numpy()
    zbar_1995 = df[f"zbar_1995_{num_sites}Sites"].to_numpy()
    forest_area_1995 = df[f"forestArea_1995_ha_{num_sites}Sites"].to_numpy()

    # Normalize Z data
    zbar_1995 /= norm_fac
    z_1995 /= norm_fac
    forest_area_1995 /= norm_fac

    # Read municipal level data
    municipal_theta_df = gpd.read_file(data_dir / "muni_data_theta.geojson")
    municipal_gamma_df = gpd.read_file(data_dir / "muni_data_gamma.geojson")

    # Read site level data
    site_theta_2017 = gpd.read_file(data_dir / f"site_{num_sites}_data_theta.geojson")
    site_gamma_2017 = gpd.read_file(data_dir / f"site_{num_sites}_data_gamma.geojson")

    # Remove geometries
    site_theta_2017_df = site_theta_2017.iloc[:, :-1]
    site_gamma_2017_df = site_gamma_2017.iloc[:, :-1]

    # Get theta and gamma samples
    baseline_fit = baseline.sample(
        num_sites=num_sites, iter_sampling=10**4, chains=5, seed=1
    )

    theta = baseline_fit.stan_variable("theta").mean(axis=0)
    gamma = baseline_fit.stan_variable("gamma").mean(axis=0)

    print(f"Data successfully loaded from {data_dir}")
    return (
        zbar_1995,
        z_1995,
        forest_area_1995,
        site_theta_2017_df,
        site_gamma_2017_df,
        municipal_theta_df,
        municipal_gamma_df,
        z_2008,
        theta,
        gamma,
    )
