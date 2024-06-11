import os
import numpy as np
import pandas as pd
import argparse
from pysrc.sampling import mpc
from pysrc.services.file_service import get_path
## This script will create markov chain samples for both constrained and unconstrained model

def mpc_simulation():
    num_site=78

    # Define the root folders
    root_folder = str(get_path("gams_file"))+"/mpc"

    computation_folder=os.path.join(root_folder, 'computation')
    simulation_folder=os.path.join(root_folder, 'simulation')
    calculation_folder_unconstrained = os.path.join(root_folder+'/computation', 'model_unconstrained')
    mc_samples_folder_unconstrained = os.path.join(root_folder+'/simulation', 'sample_unconstrained')
    calculation_folder_constrained = os.path.join(root_folder+'/computation', 'model_constrained')
    mc_samples_folder_constrained = os.path.join(root_folder+'/simulation', 'sample_constrained')


    os.makedirs(computation_folder, exist_ok=True)
    os.makedirs(simulation_folder, exist_ok=True)
    os.makedirs(calculation_folder_unconstrained, exist_ok=True)
    os.makedirs(mc_samples_folder_unconstrained, exist_ok=True)
    os.makedirs(calculation_folder_constrained, exist_ok=True)
    os.makedirs(mc_samples_folder_constrained, exist_ok=True)

    # get mc samples
    mpc.mc_samples_unconstrained(location=mc_samples_folder_unconstrained)
    mpc.mc_samples_constrained(location=mc_samples_folder_constrained)
    print("sampling is done")


    # get gdx files
    mpc.gdx_files(num_sites=78,location=calculation_folder_unconstrained)
    mpc.gdx_files(num_sites=78,location=calculation_folder_constrained)

    print("gdx file is done")



    # move samples into computation folders
    mpc.paste_file(root=root_folder,ori=mc_samples_folder_unconstrained,des=calculation_folder_unconstrained)
    mpc.paste_file(root=root_folder,ori=mc_samples_folder_constrained,des=calculation_folder_constrained,model="constrained")


    print("file movement is done")

    return print("mpc simulation is done")


