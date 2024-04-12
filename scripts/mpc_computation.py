import os
import time

import numpy as np
import pandas as pd
from gams import GamsWorkspace


import argparse
parser = argparse.ArgumentParser(description="parameter settings")
parser.add_argument("--id",type=int,default=1)
parser.add_argument("--pf",type=float,default=20.76)
parser.add_argument("--mc", type=str, default="unconstrained")


args = parser.parse_args()
pf = args.pf
id = args.id
mc_model = args.mc




# set the path
root_folder = os.getcwd()+"/gams_file/mpc"
calculation_folder_unconstrained = os.path.join(root_folder+'/computation', 'model_unconstrained')
calculation_folder_constrained = os.path.join(root_folder+'/computation', 'model_constrained')

# set the model
if mc_model == 'constrained':
    comp_path = calculation_folder_constrained
    gams_file = "mpc_con.gms"
else:
    comp_path = calculation_folder_unconstrained
    gams_file = "mpc_uncon.gms"
    
working_directory = comp_path + f"/mc_{id}/"

# Create Gams Workspace and run gms files

ws = GamsWorkspace(
    system_directory=os.getcwd()
    + "/gams_file/gams45.1_linux_x64_64_sfx",
    working_directory = comp_path + f"/mc_{id}/",
)

db = ws.add_database(in_model_name="myDB")
db.add_parameter("p_e", 0).add_record().value = pf


t1 = ws.add_job_from_file(gams_file)
t1.run(databases=db)

   