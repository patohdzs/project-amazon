import os
import pickle
import shutil

from pysrc.optimization import solve_planner_problem
from pysrc.sampling import baseline
from pysrc.services.data_service import load_site_data
from pysrc.services.file_service import get_path


def get_optimization(
    solver="gurobi",
    num_sites=78,
    pee=5,
    pa=41.11,
    model="det",
    xi=1,
):
    (
        zbar_2017,
        z_2017,
        forest_area_2017,
        _,
        _,
        _,
        _,
    ) = load_site_data(num_sites)

    if model == "det":
        # Set initial theta & gamma using baseline mean
        baseline_fit = baseline.sample(
            num_sites=num_sites, iter_sampling=10**4, chains=5, seed=1
        )
        theta_vals = baseline_fit.stan_variable("theta").mean(axis=0)
        gamma_vals = baseline_fit.stan_variable("gamma").mean(axis=0)
        x0_vals = gamma_vals * forest_area_2017

        b = [0, 10, 15, 20, 25]
        pe_values = [pee + bi for bi in b]
        for pe in pe_values:
            solve_planner_problem(
                time_horizon=200,
                theta=theta_vals,
                gamma=gamma_vals,
                x0=x0_vals,
                zbar=zbar_2017,
                z0=z_2017,
                price_emissions=pe,
                price_cattle=pa,
            )
            print("Results for pe = ", pe)
            output_base_path = str(get_path("output"))
            gams_working_directory = str(get_path("gams_file")) + f"/{num_sites}sites/"
            subfolder_path = os.path.join(
                output_base_path,
                "optimization",
                model,
                solver,
                f"{num_sites}sites",
                f"pa_{pa}",
                f"pe_{pe}",
            )

            if not os.path.exists(subfolder_path):
                os.makedirs(subfolder_path)

            file_names = [
                "amazon_data_z.dat",
                "amazon_data_x.dat",
                "amazon_data_u.dat",
                "amazon_data_v.dat",
                "amazon_data_w.dat",
            ]
            for file_name in file_names:
                source_file = os.path.join(gams_working_directory, file_name)
                destination_file = os.path.join(subfolder_path, file_name)
                shutil.copy(source_file, destination_file)

    if model == "hmc":
        result_folder = os.path.join(
            str(get_path("output")),
            "sampling",
            solver,
            f"{num_sites}sites",
            f"pa_{pa}",
            f"xi_{xi}",
        )
        b = [0, 10, 15, 20, 25]
        pe_values = [pee + bi for bi in b]
        for pe in pe_values:
            with open(result_folder + f"/pe_{pe}/results.pcl", "rb") as f:
                b = pickle.load(f)
            theta_vals = b["final_sample"][:16000, :78].mean(axis=0)
            gamma_vals = b["final_sample"][:16000, 78:].mean(axis=0)
            x0_vals = gamma_vals * forest_area_2017

            solve_planner_problem(
                time_horizon=200,
                theta=theta_vals,
                gamma=gamma_vals,
                x0=x0_vals,
                zbar=zbar_2017,
                z0=z_2017,
                price_emissions=pe,
                price_cattle=pa,
            )
            print("Results for pe = ", pe)
            output_base_path = str(get_path("output"))
            gams_working_directory = str(get_path("gams_file")) + f"/{num_sites}sites/"
            subfolder_path = os.path.join(
                output_base_path,
                "optimization",
                model,
                solver,
                f"{num_sites}sites",
                f"pa_{pa}",
                f"pe_{pe}",
            )

            if not os.path.exists(subfolder_path):
                os.makedirs(subfolder_path)

            file_names = [
                "amazon_data_z.dat",
                "amazon_data_x.dat",
                "amazon_data_u.dat",
                "amazon_data_v.dat",
                "amazon_data_w.dat",
            ]
            for file_name in file_names:
                source_file = os.path.join(gams_working_directory, file_name)
                destination_file = os.path.join(subfolder_path, file_name)
                shutil.copy(source_file, destination_file)
    return
