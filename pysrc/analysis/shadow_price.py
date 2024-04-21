import numpy as np
import pandas as pd
import geopandas as gpd
from pathlib import Path
from pysrc.optimization import gams, gurobi
from pysrc.sampling import baseline
from pysrc.services.file_service import get_path


def load_site_data_1995(
    num_sites: int,
    norm_fac: float = 1e9,
    data_dir: Path = get_path("data", "hmc"),
):
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
    
    
    baseline_fit = baseline.sample(
        model_name="full_model",
        num_sites=num_sites,
        iter_sampling=10**4,
        chains=5,
        seed=1,
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



def shadow_price_opt(
    zbar_1995,
    z_1995,
    forest_area_1995,
    z_2008,
    theta,
    gamma,
    sitenum = 78,
    opt = 'gams',
    timehzn = 200,
    pa=41.11,
    pe=7.1,
):
    
    # Computing carbon absorbed in start period
    x0_vals_1995 = gamma * forest_area_1995

    # Choose optimizer
    if opt == "gurobi":
        solve_planner_problem = gurobi.solve_planner_problem

    elif opt == "gams":
        solve_planner_problem = gams.solve_planner_problem

    else:
        raise ValueError("Optimizer must be one of ['gurobi', 'gams']")

    results = solve_planner_problem(
        T=timehzn,
        theta=theta,
        gamma=gamma,
        x0=x0_vals_1995,
        zbar=zbar_1995,
        z0=z_1995,
        pe=pe,
        pa=pa,
        model="shadow_price",
    )
    Z=results['Z']
    z_2008_agg=np.sum(z_2008)/1e9
    ratio=np.abs((np.sum(Z[13])-z_2008_agg)/z_2008_agg)
    
    return ratio


def shadow_price_cal(sitenum = 78,
                     pa=41.11,
                     opt='gams',):
    (
        zbar_1995,
        z_1995,
        forest_area_1995,
        _,
        _,
        _,
        _,
        z_2008,
        theta,
        gamma
    ) = load_site_data_1995(sitenum)
        
    
    pe_values = np.arange(5, 10, 0.1)
    results = np.array([shadow_price_opt(    zbar_1995,
                                         z_1995,
                                         forest_area_1995,
                                         z_2008,
                                         theta,
                                         gamma,
                                         sitenum=sitenum,
                                         opt=opt,
                                         timehzn=200,
                                         pa=pa,
                                         pe=pe,
                                         ) for pe in pe_values])
    
    
    min_index = np.argmin(results)
    min_result = results[min_index]
    min_pe = pe_values[min_index]
    

    return min_result,min_pe

    

    
min_result,min_pe=shadow_price_cal(sitenum=1043)
print("min_result",min_result,"min_pe",min_pe)