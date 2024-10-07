import os
import pickle
import shutil
import numpy as np
import pandas as pd

from pysrc.optimization import solve_planner_problem
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
    ) = load_site_data(num_sites)

    if model == "det":
        # Set initial theta & gamma using baseline mean
        # baseline_fit = baseline.sample(
        #     num_sites=num_sites, iter_sampling=10**4, chains=5, seed=1
        # )
        theta_vals = (
            pd.read_csv(
                get_path("data", "calibration", "hmc") / f"theta_fit_{num_sites}.csv"
            )
            .to_numpy()[:,]
            .flatten()
        )
        gamma_vals = (
            pd.read_csv(
                get_path("data", "calibration", "hmc") / f"gamma_fit_{num_sites}.csv"
            )
            .to_numpy()[:,]
            .flatten()
        )
        x0_vals = gamma_vals * forest_area_2017

        b = [0, 10, 15, 20, 25]
        pe_values = [pee + bi for bi in b]
        for pe in pe_values:
            results=solve_planner_problem(
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

            output_folder = os.path.join(
                str(get_path("output")),
                "optimization",
                model,
                solver,
                f"{num_sites}sites",
                f"pa_{pa}",
                f"pe_{pe}",
            )

            if not os.path.exists(output_folder):
                os.makedirs(output_folder)
                
            np.savetxt(os.path.join(output_folder, "Z.txt"), results.Z, delimiter=",")
            np.savetxt(os.path.join(output_folder, "X.txt"), results.X, delimiter=",")
            np.savetxt(os.path.join(output_folder, "U.txt"), results.U, delimiter=",")
            np.savetxt(os.path.join(output_folder, "V.txt"), results.V, delimiter=",")


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

            results=solve_planner_problem(
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

            output_folder = os.path.join(
                str(get_path("output")),
                "optimization",
                model,
                solver,
                f"{num_sites}sites",
                f"pa_{pa}",
                f"pe_{pe}",
            )

            if not os.path.exists(output_folder):
                os.makedirs(output_folder)
                
            np.savetxt(os.path.join(output_folder, "Z.txt"), results.Z, delimiter=",")
            np.savetxt(os.path.join(output_folder, "X.txt"), results.X, delimiter=",")
            np.savetxt(os.path.join(output_folder, "U.txt"), results.U, delimiter=",")
            np.savetxt(os.path.join(output_folder, "V.txt"), results.V, delimiter=",")
    return
