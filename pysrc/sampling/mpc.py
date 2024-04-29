import os
import shutil

import numpy as np
import pandas as pd

from pysrc.sampling import baseline, hmm_estimation

from ..services.data_service import load_site_data


def sample_price_paths(location, var="uncon"):
    np.random.seed(123)

    price_states, M = hmm_estimation.estimate_price_model(num_iter=5, var=var)

    num_simulations = 200

    states = {"low": price_states[0], "high": price_states[1]}
    initial_state = "high"

    probability_matrix = {
        "low": {"low": M[0, 0], "high": M[0, 1]},
        "high": {"low": M[1, 0], "high": M[1, 1]},
    }

    # Number of observations
    num_observations = 200

    for i in range(1, num_simulations + 1):
        # Generating the Markov chain for each simulation
        current_state = initial_state
        markov_chain = [states[current_state]]

        for _ in range(num_observations - 1):
            next_state = np.random.choice(
                list(probability_matrix[current_state].keys()),
                p=list(probability_matrix[current_state].values()),
            )
            markov_chain.append(states[next_state])
            current_state = next_state

        transformed_markov_chain = [
            1 if price == price_states[0] else 2 for price in markov_chain
        ]

        # Creating an index from 1 to 200
        index = range(1, num_observations + 1)

        # Combine index and transformed Markov chain into a DataFrame
        markov_chain_df = pd.DataFrame(
            {"Index": index, "scenario": transformed_markov_chain}
        )

        # Specify the filename for each CSV file
        csv_filename = f"/mc_{i}.csv"
        output = location
        # Output to CSV file
        markov_chain_df.to_csv(output + csv_filename, index=False)

    return "mc sampling is done"


def gdx_files(location, num_sites=78):
    # Load site data

    (
        zbar_2017,
        z_2017,
        forest_area_2017,
        _,
        _,
        _,
        _,
    ) = load_site_data(num_sites)

    baseline_fit = baseline.sample(
        model_name="full_model",
        num_sites=num_sites,
        iter_sampling=10**4,
        chains=5,
        seed=1,
    )

    theta = baseline_fit.stan_variable("theta").mean(axis=0)
    gamma = baseline_fit.stan_variable("gamma").mean(axis=0)

    scale = 1e9
    # Computing carbon absorbed in start period
    x0_vals = gamma * forest_area_2017

    x0_vals = x0_vals * scale
    site_z_vals = z_2017 * scale
    zbar_2017 = zbar_2017 * scale

    x0data = pd.DataFrame(x0_vals)
    x0data.columns = ["x_2017_78Sites"]
    x0data.index = x0data.index + 1
    saveto = os.path.join(location, "X0Data.csv")
    x0data.to_csv(saveto)

    z0data = pd.DataFrame(site_z_vals)
    z0data.columns = ["z_2017_78Sites"]
    z0data.index = z0data.index + 1
    saveto = os.path.join(location, "Z0Data.csv")
    z0data.to_csv(saveto)

    zbardata = pd.DataFrame(zbar_2017)
    zbardata.columns = ["zbar_2017_78Sites"]
    zbardata.index = zbardata.index + 1
    saveto = os.path.join(location, "ZbarData.csv")
    zbardata.to_csv(saveto)

    gammadata = pd.DataFrame(gamma)
    gammadata.columns = ["gamma_78Sites"]
    gammadata.index = gammadata.index + 1
    saveto = os.path.join(location, "GammaData.csv")
    gammadata.to_csv(saveto)

    thetadata = pd.DataFrame(theta)
    thetadata.columns = ["theta_78Sites"]
    thetadata.index = thetadata.index + 1
    saveto = os.path.join(location, "ThetaData.csv")
    thetadata.to_csv(saveto)

    return "gdx file is done"


def paste_file(root, ori, des, model="unconstrained"):
    # Create subfolders 'mc_1' to 'mc_100' inside the 'calculation' folder
    for i in range(1, 201):
        mc_subfolder = os.path.join(des, f"mc_{i}")
        os.makedirs(mc_subfolder, exist_ok=True)

    # Copy 'mc_i.csv' files into their respective subfolders
    for i in range(1, 201):
        mc_csv_file_src = os.path.join(ori, f"mc_{i}.csv")
        mc_csv_file_dest = os.path.join(des, f"mc_{i}", "mc_1.csv")
        shutil.copy(mc_csv_file_src, mc_csv_file_dest)

    # Copy the 'gms' file into each subfolder
    if model == "constrained":
        target_file = "mpc_con.gms"
    else:
        target_file = "mpc_uncon.gms"

    gms_file = os.path.join(root, target_file)
    for i in range(1, 201):
        mc_subfolder = os.path.join(des, f"mc_{i}")
        shutil.copy(gms_file, mc_subfolder)

    gms_file = os.path.join(des, "GammaData.csv")
    for i in range(1, 201):
        mc_subfolder = os.path.join(des, f"mc_{i}")
        shutil.copy(gms_file, mc_subfolder)

    gms_file = os.path.join(des, "ThetaData.csv")
    for i in range(1, 201):
        mc_subfolder = os.path.join(des, f"mc_{i}")
        shutil.copy(gms_file, mc_subfolder)

    gms_file = os.path.join(des, "X0Data.csv")
    for i in range(1, 201):
        mc_subfolder = os.path.join(des, f"mc_{i}")
        shutil.copy(gms_file, mc_subfolder)

    gms_file = os.path.join(des, "Z0Data.csv")
    for i in range(1, 201):
        mc_subfolder = os.path.join(des, f"mc_{i}")
        shutil.copy(gms_file, mc_subfolder)

    gms_file = os.path.join(des, "ZbarData.csv")
    for i in range(1, 201):
        mc_subfolder = os.path.join(des, f"mc_{i}")
        shutil.copy(gms_file, mc_subfolder)

    return "movement is done"
