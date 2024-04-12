import os
import time

import numpy as np
import pandas as pd
from gams import GamsWorkspace
from ..services.file_service import get_path



def solve_planner_problem(
    T,
    theta,
    gamma,
    x0,
    zbar,
    z0,
    alpha=0.045007414,
    kappa=2.094215255,
    pe=20.76,
    pa=44.75,
    zeta=1.66e-4 * 1e9,  # use the same normalization factor
    scale=1e9,
):
    # Default Path to the data folder
    _DATA_DIR = str(get_path("gams_file"))
    num_sites = theta.shape[0]

    x0_vals = x0 * scale
    site_z_vals = z0 * scale
    zbar_2017 = zbar * scale

    working_directory = _DATA_DIR + f"/{num_sites}sites/"

    x0data = pd.DataFrame(x0_vals)
    x0data.index = x0data.index + 1
    saveto = os.path.join(working_directory, "X0Data.csv")
    x0data.to_csv(saveto)

    z0data = pd.DataFrame(site_z_vals)
    z0data.index = z0data.index + 1
    saveto = os.path.join(working_directory, "Z0Data.csv")
    z0data.to_csv(saveto)

    zbardata = pd.DataFrame(zbar_2017)
    zbardata.index = zbardata.index + 1
    saveto = os.path.join(working_directory, "ZbarData.csv")
    zbardata.to_csv(saveto)

    gammadata = pd.DataFrame(gamma)
    gammadata.index = gammadata.index + 1
    saveto = os.path.join(working_directory, "GammaData.csv")
    gammadata.to_csv(saveto)

    thetadata = pd.DataFrame(theta)
    thetadata.index = thetadata.index + 1
    saveto = os.path.join(working_directory, "ThetaData.csv")
    thetadata.to_csv(saveto)

    # Create Gams Workspace

    ws = GamsWorkspace(
        system_directory=_DATA_DIR
        + "/gams45.1_linux_x64_64_sfx",
        working_directory=_DATA_DIR + f"/{num_sites}sites/",
    )

    db = ws.add_database(in_model_name="myDB")
    db.add_parameter("p_e", 0).add_record().value = pe

    start_time = time.time()
    gams_file = f"hmc_{num_sites}sites.gms"
    # shutil.copy(gams_file, working_directory)
    t1 = ws.add_job_from_file(gams_file)
    t1.run(databases=db)

    readfrom = os.path.join(working_directory, "amazon_data_u.dat")
    U = pd.read_csv(readfrom, delimiter="\t").drop("T/R ", axis=1).to_numpy()


    readfrom = os.path.join(working_directory, "amazon_data_v.dat")
    V = pd.read_csv(readfrom, delimiter="\t").drop("T/R ", axis=1).to_numpy()


    readfrom = os.path.join(working_directory, "amazon_data_w.dat")
    dfw = pd.read_csv(readfrom, delimiter="\t")
    dfw = dfw.drop("T   ", axis=1)
    w = dfw.to_numpy()[:,0]
    # os.remove(readfrom)

    readfrom = os.path.join(working_directory, "amazon_data_x.dat")
    dfx = pd.read_csv(readfrom, delimiter="\t")
    dfx = dfx.drop("T   ", axis=1)
    X = dfx.to_numpy()


    readfrom = os.path.join(working_directory, "amazon_data_z.dat")
    dfz = pd.read_csv(readfrom, delimiter="\t").drop("T/R ", axis=1)
    Z = dfz.to_numpy()


    print(f"Done! Time elapsed: {time.time()-start_time} seconds.")

    return {
        "Z": Z,
        "X": X,
        "U": U,
        "V": V,
        "w": w,
    }



def vectorize_trajectories(Z, X, U, V, w):
    X_agg = X.sum(axis=1)
    X_agg = X_agg.reshape(X_agg.size, 1)

    sol_val_Ua = (w[:-1] ** 2).T.flatten()
    sol_val_X = np.concatenate((Z.T, X_agg.T, np.ones((1, Z.T.shape[1]))))
    sol_val_Up = U[:-1, :].T
    sol_val_Um = V[:-1, :].T
    sol_val_Z = sol_val_Up - sol_val_Um
    return {
        "sol_val_X": sol_val_X,
        "sol_val_Up": sol_val_Up,
        "sol_val_Um": sol_val_Um,
        "sol_val_Z": sol_val_Z,
        "sol_val_Ua": sol_val_Ua,
    }