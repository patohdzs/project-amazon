import os
import numpy as np
import pandas as pd
import argparse
from pysrc.sampling import mpc


# Read arguments from stdin
parser = argparse.ArgumentParser(description="parameter settings")
parser.add_argument("--sitenum", type=int, default=78)
# Parse arguments
args = parser.parse_args()


# Define the root folders
root_folder = os.getcwd()+"/gams_file/mpc"
calculation_folder_unconstrained = os.path.join(root_folder, 'model_unconstrained')
mc_samples_folder_unconstrained = os.path.join(root_folder, 'sample_unconstrained')
calculation_folder_constrained = os.path.join(root_folder, 'model_constrained')
mc_samples_folder_constrained = os.path.join(root_folder, 'sample_constrained')


# get mc samples
mpc.mc_samples_unconstrained(location=mc_samples_folder_unconstrained)
mpc.mc_samples_constrained(location=mc_samples_folder_constrained)
print("sampling is done")


# get gdx files
mpc.gdx_files(num_sites=args.sitenum,location=calculation_folder_unconstrained)
mpc.gdx_files(num_sites=args.sitenum,location=calculation_folder_constrained)

print("gdx file is done")





# # Create subfolders 'mc_1' to 'mc_100' inside the 'calculation' folder
# for i in range(1, 201):
#     mc_subfolder = os.path.join(calculation_folder_unconstrained, f'mc_{i}')
#     os.makedirs(mc_subfolder, exist_ok=True)

# # Copy 'mc_i.csv' files into their respective subfolders
# for i in range(1, 201):
#     mc_csv_file_src = os.path.join(mc_samples_folder_unconstrained, f'mc_{i}.csv')
#     mc_csv_file_dest = os.path.join(calculation_folder_unconstrained, f'mc_{i}', 'mc_1.csv')
#     shutil.copy(mc_csv_file_src, mc_csv_file_dest)

# # Copy the 'gms' file into each subfolder
# gms_file = os.path.join(calculation_folder_unconstrained, 'mpc_78sites_4_model3.gms')
# for i in range(1, 201):
#     mc_subfolder = os.path.join(calculation_folder_unconstrained, f'mc_{i}')
#     shutil.copy(gms_file, mc_subfolder)


# gms_file = os.path.join(calculation_folder_unconstrained, 'GammaData.csv')
# for i in range(1, 201):
#     mc_subfolder = os.path.join(calculation_folder_unconstrained, f'mc_{i}')
#     shutil.copy(gms_file, mc_subfolder)
    
# gms_file = os.path.join(calculation_folder_unconstrained, 'ThetaData.csv')
# for i in range(1, 201):
#     mc_subfolder = os.path.join(calculation_folder_unconstrained, f'mc_{i}')
#     shutil.copy(gms_file, mc_subfolder)
    
    
# gms_file = os.path.join(calculation_folder_unconstrained, 'X0Data.csv')
# for i in range(1, 201):
#     mc_subfolder = os.path.join(calculation_folder_unconstrained, f'mc_{i}')
#     shutil.copy(gms_file, mc_subfolder)
    
# gms_file = os.path.join(calculation_folder_unconstrained, 'Z0Data.csv')
# for i in range(1, 201):
#     mc_subfolder = os.path.join(calculation_folder_unconstrained, f'mc_{i}')
#     shutil.copy(gms_file, mc_subfolder)
    
    
# gms_file = os.path.join(calculation_folder_unconstrained, 'ZbarData.csv')
# for i in range(1, 201):
#     mc_subfolder = os.path.join(calculation_folder_unconstrained, f'mc_{i}')
#     shutil.copy(gms_file, mc_subfolder)




