import os
import time
import shutil
import numpy as np
import pandas as pd
from gams import GamsWorkspace
from ..services.file_service import get_path
from ..services.data_service import load_price_data


def solve_planner_problem(
    T,
    theta,
    gamma,
    x0,
    zbar,
    z0,
    alpha=0.045007414,
    delta=0.02,
    kappa=2.094215255,
    dt=1,
    pe=20.76,
    pa=44.75,
    zeta=1.66e-4 * 1e9,  # use the same normalization factor
    scale=1e9,
    model="prediction",
):
    # Default Path to the data folder
    _DATA_DIR = str(get_path("gams_file"))
    num_sites = theta.shape[0]

    x0_vals = x0 * scale
    site_z_vals = z0 * scale
    zbar_2017 = zbar * scale

    working_directory = _DATA_DIR + f"/{num_sites}sites/"

    if model=="prediction":
        pa_data = pd.DataFrame({
        'price': [pa] * 200  
        })
        pa_data.index = range(1, 201)
        pa_data.index.name = None
        saveto = os.path.join(working_directory, "p_a.csv")
        pa_data.to_csv(saveto, header=True, index_label=None)
    elif model=="shadow_price":
        pa_list=load_price_data()
        prices = np.concatenate((pa_list, np.full(200 - len(pa_list), pa)))
        pa_data = pd.DataFrame({
            'price': prices
        })
        pa_data.index = range(1, 201)
        pa_data.index.name = None
        saveto = os.path.join(working_directory, "p_a.csv")
        pa_data.to_csv(saveto, header=True, index_label=None)

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

    # ws = GamsWorkspace(
    #     system_directory=os.path.dirname(os.path.dirname(os.getcwd()))
    #     + "/gams/gams45.1_linux_x64_64_sfx",
    #     working_directory=os.getcwd() + f"/gams_file/{num_sites}sites/",
    # )
    
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
    
    
    
    
    




def mpc_shadow_price(
    T,
    theta,
    gamma,
    x0,
    zbar,
    z0,
    alpha=0.045007414,
    delta=0.02,
    kappa=2.094215255,
    dt=1,
    pe=20.76,
    pa=44.75,
    zeta=1.66e-4 * 1e9,  # use the same normalization factor
    scale=1e9,
    model = 'shadow_price',
):
    # Default Path to the data folder
    _DATA_DIR = str(get_path("gams_file"))
    num_sites = theta.shape[0]

    x0_vals = x0 * scale
    site_z_vals = z0 * scale
    zbar_2017 = zbar * scale

    working_directory = _DATA_DIR + "/mpc/shadowprice/"
    if not os.path.exists(working_directory):
        os.makedirs(working_directory)

    gms_file = _DATA_DIR+"/mpc/mpc_shadowprice.gms"
    shutil.copy(gms_file, working_directory)
    gams_file = "mpc_shadowprice.gms"
    

        
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
    
    ws = GamsWorkspace(

        system_directory=_DATA_DIR
        + "/gams45.1_linux_x64_64_sfx",
        working_directory=working_directory,

    )

    db = ws.add_database(in_model_name="myDB")
    db.add_parameter("p_e", 0).add_record().value = pe

    start_time = time.time()
    t1 = ws.add_job_from_file(gams_file)
    t1.run(databases=db)
    
    data = read_dat_file(working_directory+'amazon_data_z.dat')

    Z = np.array([np.array(data[i]['values']) for i in range(1, 20)])/1e9


    print(f"Done! Time elapsed: {time.time()-start_time} seconds.")

    return {
        "Z": Z,
    }
    
    
    
    
    
def read_dat_file(filename):
    # Create an empty dictionary to store our data.
    data_dict = {}
    
    # Open and read the file.
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Initialize variables to temporarily store the current Y value and the list of numbers.
    current_y = None
    values = []
    p_a = None

    # Iterate over each line in the file.
    for line in lines:
        line = line.strip()  # Remove whitespace.

        # If the line starts with 'Y =', we've reached a new block of data.
        if line.startswith('Y ='):
            # If we have a previous block of data, store it in our dictionary.
            if current_y is not None:
                data_dict[current_y] = {'p_a': p_a, 'values': values}
                values = []

            # Split and store the Y and p_a values.
            parts = line.split()
            current_y = int(parts[2])
            # p_a = float(parts[5])

        # Otherwise, assume the line contains a number and add it to our list of values.
        else:
            try:
                value = float(line)
                values.append(value)
            except ValueError:
                pass

    # Add the last block of data to our dictionary.
    if current_y is not None:
        data_dict[current_y] = {'p_a': p_a, 'values': values}

    return data_dict