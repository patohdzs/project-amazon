import os
import shutil
import time
from concurrent.futures import ProcessPoolExecutor
from gams import GamsWorkspace
from pysrc.services.file_service import get_path

# Function to process each ID
def process_id(id, pf=16.9, mc_model='unconstrained'):
    root_folder = os.getcwd() + "/gams_file/mpc"
    comp_path = os.path.join(root_folder, 'computation', 'model_constrained' if mc_model == 'constrained' else 'model_unconstrained')
    gams_file = "mpc_con.gms" if mc_model == 'constrained' else "mpc_uncon.gms"
    
    # Setup working directory for each ID
    working_directory = os.path.join(comp_path, f"mc_{id}")

    # Create GAMS workspace and run GAMS file
    ws = GamsWorkspace(
        system_directory=os.getcwd() + "/gams_file/gams45.1_linux_x64_64_sfx", working_directory=working_directory)
    db = ws.add_database(in_model_name="myDB")
    db.add_parameter("p_e", 0).add_record().value = pf
    t1 = ws.add_job_from_file(gams_file)
    t1.run(databases=db)

    
    # Copy output files
    output_base_path = get_path("output")
    gams_working_directory = working_directory
    subfolder_path = os.path.join(output_base_path, "optimization", "mpc", "gams", '78sites', f'model_{mc_model}', f'pe_{pf}', f'mc_{id}')
    if not os.path.exists(subfolder_path):
        os.makedirs(subfolder_path)
    file_names = ['amazon_data_z.dat', 'amazon_data_x.dat', 'amazon_data_u.dat', 'amazon_data_v.dat', 'amazon_data_w.dat']
    for file_name in file_names:
        source_file = os.path.join(gams_working_directory, file_name)
        destination_file = os.path.join(subfolder_path, file_name)
        shutil.copy(source_file, destination_file)

# Example usage of concurrent processing

def parallel_conduction(idstart=1,idend=5,pf=6.9,model="unconstrained"):
    ids = range(idstart,idend+1)  
    start_time = time.time()

    with ProcessPoolExecutor(max_workers=5) as executor:
        executor.map(process_id, ids, [pf]*len(ids), [model]*len(ids))
        
    total_time = time.time() - start_time
    return print(f"Id{idstart} to Id{idend} done! Execution time: {total_time:.2f} seconds")
