import os
import time

import numpy as np
import pandas as pd
from gams import GamsWorkspace

import shutil

from pysrc.services.file_service import get_path


import argparse
parser = argparse.ArgumentParser(description="parameter settings")
parser.add_argument("--id",type=int,default=1)
parser.add_argument("--pe",type=float,default=20.76)


args = parser.parse_args()

pe = args.pe
id = args.id



workdir = os.getcwd()

print("mc_file",id)    
mc_model='unconstrained'
    

root_folder = os.getcwd() + "/gams_file/mpc"
comp_path = os.path.join(root_folder, 'computation', 'model_constrained' if mc_model == 'constrained' else 'model_unconstrained')
gams_file = "mpc_con.gms" if mc_model == 'constrained' else "mpc_uncon.gms"

# Setup working directory for each ID
working_directory = os.path.join(comp_path, f"mc_{id}")

# Create GAMS workspace and run GAMS file
ws = GamsWorkspace(
    system_directory=os.getcwd() + "/gams_file/gams45.1_linux_x64_64_sfx", working_directory=working_directory)
db = ws.add_database(in_model_name="myDB")
db.add_parameter("p_e", 0).add_record().value = pe
t1 = ws.add_job_from_file(gams_file)
t1.run(databases=db)


# Copy output files
output_base_path = get_path("output")
gams_working_directory = working_directory
subfolder_path = os.path.join(output_base_path, "optimization", "mpc", "gams", '78sites', f'model_{mc_model}', f'pe_{pe}', f'mc_{id}')
if not os.path.exists(subfolder_path):
    os.makedirs(subfolder_path)
file_names = ['amazon_data_z.dat', 'amazon_data_x.dat', 'amazon_data_u.dat', 'amazon_data_v.dat', 'amazon_data_w.dat']
for file_name in file_names:
    source_file = os.path.join(gams_working_directory, file_name)
    destination_file = os.path.join(subfolder_path, file_name)
    shutil.copy(source_file, destination_file)
    
    
print("all done")