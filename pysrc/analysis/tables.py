import os
from pysrc.services.file_service import get_path
import numpy as np
import pandas as pd
from pysrc.sampling import baseline
import pickle

def format_float(value):
    return f"{value:.2f}"

def read_theta(num_sites):
    baseline_fit = baseline.sample(
    model_name="full_model",
    num_sites=num_sites,
    iter_sampling=10**4,
    chains=5,
    seed=1,
    )
        
    dft_np = baseline_fit.stan_variable("theta").mean(axis=0)
    return dft_np

def read_file(result_directory):
    dfz = pd.read_csv(result_directory + "/amazon_data_z.dat", delimiter="\t")
    dfz = dfz.drop("T/R ", axis=1)
    dfz_np = dfz.to_numpy()

    dfx = pd.read_csv(result_directory + "/amazon_data_x.dat", delimiter="\t")
    dfx = dfx.drop("T   ", axis=1)
    dfx_np = dfx.to_numpy()
    dfxdot = np.diff(dfx_np, axis=0)

    dfu = pd.read_csv(result_directory + "/amazon_data_u.dat", delimiter="\t")
    dfu = dfu.drop("T/R ", axis=1)
    dfu_np = dfu.to_numpy()

    dfv = pd.read_csv(result_directory + "/amazon_data_v.dat", delimiter="\t")
    dfv = dfv.drop("T/R ", axis=1)
    dfv_np = dfv.to_numpy()

    return (dfz_np/1e2, dfxdot/1e2, dfu_np/1e2, dfv_np/1e2)

def value_decom(pee=7.1,num_sites=78,opt='gams',pa=41.11,model='det',xi=1):
    
    b = [0, 10, 15, 20, 25]
    pe = [pee + bi for bi in b]
    kappa = 2.094215255
    zeta = 1.66e-4 * 1e11
    output_folder = str(get_path("output"))+"/tables/"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if model == 'det':
        dft_np = read_theta(num_sites)

        
    results = []
    for j in range(5):
        order = j

        result_folder =  os.path.join(str(get_path("output")), "optimization",model,opt,f'{num_sites}sites',f'pa_{pa}' ,f'pe_{pe[order]}')

        if model == 'hmc':
            theta_folder = os.path.join(str(get_path("output")), "sampling",opt,f'{num_sites}sites',f'pa_{pa}' ,f'xi_{xi}')
            with open(theta_folder+f"/pe_{pe[order]}/results.pcl", 'rb') as f:
                para_file = pickle.load(f)
            dft_np=para_file['final_sample'][:16000,:78].mean(axis=0)
            
        (dfz_np, dfxdot, dfu_np, dfv_np) = read_file(result_folder)

        results_AO = []
        for i in range(200):
            result_AO = pa * np.dot(dfz_np[i + 1], dft_np) / ((1 + 0.02) ** (i))
            results_AO.append(result_AO)
        total_AO = np.sum(results_AO)

        results_NT = []
        for i in range(200):
            result_NT = (
                -b[order]
                * (kappa * np.sum(dfz_np[i + 1]) - dfxdot[i])
                / ((1 + 0.02) ** (i))
            )
            results_NT.append(result_NT)
        total_NT = np.sum(results_NT)

        results_CS = []
        for i in range(200):
            result_CS = (
                -pee * (kappa * np.sum(dfz_np[i + 1]) - dfxdot[i]) / ((1 + 0.02) ** (i))
            )
            results_CS.append(result_CS)
        total_CS = np.sum(results_CS)

        results_AC = []
        for i in range(200):
            result_AC = (
                (zeta / 2)
                * (np.sum(dfu_np[i]) + np.sum(dfv_np[i])) ** 2
                / ((1 + 0.02) ** (i))
            )
            results_AC.append(result_AC)
        total_AC = np.sum(results_AC)

        total_PV = total_AO + total_NT + total_CS - total_AC

        total_AO = total_AO.round(2)
        total_NT = total_NT.round(2)
        total_CS = total_CS.round(2)
        total_AC = total_AC.round(2)
        total_PV = total_PV.round(2)

        iteration_results = {
            "pa": format_float(pa),
            "pe": format_float(pe[order]),
            "b": b[order],
            "agricultural output value": format_float(total_AO),
            "net transfers": format_float(total_NT),
            "forest services": format_float(total_CS),
            "adjustment costs": format_float(total_AC),
            "planner value": format_float(total_PV),
        }

        results.append(iteration_results)

    results_df = pd.DataFrame(results)
    latex_code = results_df.to_latex(index=False)


    with open(output_folder + f"present_value_site{num_sites}_pa{pa}_{model}.tex", "w") as file:
        file.write(latex_code)

    return 





def transfer_cost(pee=7.1,num_sites=78,opt='gams',pa=41.11,y=30,model='det'):
    b = [0, 10, 15, 20, 25]
    pe = [pee + bi for bi in b]
    
    baseline_folder = os.path.join(str(get_path("output")), "optimization",model,opt,f'{num_sites}sites',f'pa_{pa}' ,f'pe_{pe[0]}')
    kappa = 2.094215255
    output_folder = str(get_path("output"))+"/tables/"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)


    dfz = pd.read_csv(baseline_folder + "/amazon_data_z.dat", delimiter="\t")
    dfz = dfz.drop("T/R ", axis=1)
    dfz_np = dfz.to_numpy()/1e2

    dfx = pd.read_csv(baseline_folder + "/amazon_data_x.dat", delimiter="\t")
    dfx = dfx.drop("T   ", axis=1)
    dfx_np = dfx.to_numpy()
    dfxdot = np.diff(dfx_np, axis=0)/1e2

    results_NCE_base = []
    for i in range(y):
        result_NCE_base = -kappa * np.sum(dfz_np[i + 1]) + dfxdot[i]
        results_NCE_base.append(result_NCE_base)
    total_NCE_base = np.sum(results_NCE_base) * 100


    results_table = []
    for j in range(5):
        order = j
        
        result_folder =  os.path.join(str(get_path("output")), "optimization",model,opt,f'{num_sites}sites',f'pa_{pa}' ,f'pe_{pe[order]}')


        (dfz_np, dfxdot, dfu_np, dfv_np) = read_file(result_folder)

        results_NCE = []
        for i in range(y):
            result_NCE = -kappa * np.sum(dfz_np[i + 1]) + dfxdot[i]
            results_NCE.append(result_NCE)
        total_NCE = np.sum(results_NCE) * 100

        results_NT2 = []
        for i in range(y):
            result_NT2 = (
                -b[order]
                * (kappa * np.sum(dfz_np[i + 1]) - dfxdot[i])
                / ((1 + 0.02) ** (i))
            )

            results_NT2.append(result_NT2)
        total_NT2 = np.sum(results_NT2)

        total_EC = total_NT2 / (total_NCE - total_NCE_base) * 100

        total_NCE = total_NCE.round(2)
        total_NT2 = total_NT2.round(2)
        total_EC = total_EC.round(2)
        iteration_results = {
            "pe": format_float(pe[order]),
            "b": b[order],
            "net captured emissions": format_float(total_NCE),
            "discounted net transfers": format_float(total_NT2),
            "discounted effective costs": format_float(total_EC),
        }

        results_table.append(iteration_results)

    results_df = pd.DataFrame(results_table)
    latex_code = results_df.to_latex(index=False)

    with open(output_folder + f"transfer_cost_{num_sites}site_{pa}pa_{y}year_{model}.tex", "w") as file:
        file.write(latex_code)
    return





def ambiguity_decom(pe_det=7.1,pe_hmc=5.3,num_sites=78,opt='gams',pa=41.11,xi=1):
    
    
    kappa = 2.094215255
    zeta = 1.66e-4 * 1e11
    output_folder = str(get_path("output"))+"/tables/"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)


        
    results_det = []
    pee = pe_det
    b = [0, 10, 15, 20, 25]
    pe = [pee + bi for bi in b]
    dft_np = read_theta(num_sites)
    for j in range(5):
        order = j
        result_folder =  os.path.join(str(get_path("output")), "optimization",'det',opt,f'{num_sites}sites',f'pa_{pa}' ,f'pe_{pe[order]}')


            
        (dfz_np, dfxdot, dfu_np, dfv_np) = read_file(result_folder)

        results_AO = []
        for i in range(200):
            result_AO = pa * np.dot(dfz_np[i + 1], dft_np) / ((1 + 0.02) ** (i))
            results_AO.append(result_AO)
        total_AO = np.sum(results_AO)

        results_NT = []
        for i in range(200):
            result_NT = (
                -b[order]
                * (kappa * np.sum(dfz_np[i + 1]) - dfxdot[i])
                / ((1 + 0.02) ** (i))
            )
            results_NT.append(result_NT)
        total_NT = np.sum(results_NT)

        results_CS = []
        for i in range(200):
            result_CS = (
                -pee * (kappa * np.sum(dfz_np[i + 1]) - dfxdot[i]) / ((1 + 0.02) ** (i))
            )
            results_CS.append(result_CS)
        total_CS = np.sum(results_CS)

        results_AC = []
        for i in range(200):
            result_AC = (
                (zeta / 2)
                * (np.sum(dfu_np[i]) + np.sum(dfv_np[i])) ** 2
                / ((1 + 0.02) ** (i))
            )
            results_AC.append(result_AC)
        total_AC = np.sum(results_AC)

        total_PV = total_AO + total_NT + total_CS - total_AC



        results_det.append({
            "b": round(b[order]),
            "total_AO": total_AO,
            "total_PV": total_PV
        })
            
        
        
    results_det_df = pd.DataFrame(results_det)    



    results_hmc = []
    pee = pe_hmc
    b = [0, 10, 15, 20, 25]
    pe = [pee + bi for bi in b]
    for j in range(5):
        order = j

        result_folder =  os.path.join(str(get_path("output")), "optimization",'hmc',opt,f'{num_sites}sites',f'pa_{pa}' ,f'pe_{pe[order]}')

        theta_folder = os.path.join(str(get_path("output")), "sampling",opt,f'{num_sites}sites',f'pa_{pa}' ,f'xi_{xi}')
        with open(theta_folder+f"/pe_{pe[order]}/results.pcl", 'rb') as f:
            para_file = pickle.load(f)
        dft_np=para_file['final_sample'][:16000,:78].mean(axis=0)
            
        (dfz_np, dfxdot, dfu_np, dfv_np) = read_file(result_folder)

        results_AO = []
        for i in range(200):
            result_AO = pa * np.dot(dfz_np[i + 1], dft_np) / ((1 + 0.02) ** (i))
            results_AO.append(result_AO)
        total_AO = np.sum(results_AO)

        results_NT = []
        for i in range(200):
            result_NT = (
                -b[order]
                * (kappa * np.sum(dfz_np[i + 1]) - dfxdot[i])
                / ((1 + 0.02) ** (i))
            )
            results_NT.append(result_NT)
        total_NT = np.sum(results_NT)

        results_CS = []
        for i in range(200):
            result_CS = (
                -pee * (kappa * np.sum(dfz_np[i + 1]) - dfxdot[i]) / ((1 + 0.02) ** (i))
            )
            results_CS.append(result_CS)
        total_CS = np.sum(results_CS)

        results_AC = []
        for i in range(200):
            result_AC = (
                (zeta / 2)
                * (np.sum(dfu_np[i]) + np.sum(dfv_np[i])) ** 2
                / ((1 + 0.02) ** (i))
            )
            results_AC.append(result_AC)
        total_AC = np.sum(results_AC)

        total_PV = total_AO + total_NT + total_CS - total_AC

        results_hmc.append({
            "b": round(b[order]),
            "total_AO": total_AO,
            "total_PV": total_PV
        })
            
        
        
    results_hmc_df = pd.DataFrame(results_hmc)    


    combined_df = pd.DataFrame({
        "b": results_det_df['b'],
        "AO_det": results_det_df['total_AO'],
        "AO_hmc": results_hmc_df['total_AO'],
        "PV_det": results_det_df['total_PV'],
        "PV_hmc": results_hmc_df['total_PV']
    })

    # Calculate the percentage change for AO and PV
    combined_df['% Change AO'] = ((combined_df['AO_hmc'] - combined_df['AO_det']) / combined_df['AO_det']) * 100
    combined_df['% Change PV'] = ((combined_df['PV_hmc'] - combined_df['PV_det']) / combined_df['PV_det']) * 100

    # Display the final DataFrame
    print(combined_df)



    latex_code = combined_df.to_latex(index=False)


    with open(output_folder + f"present_value_site_ambiguity_comparison.tex", "w") as file:
        file.write(latex_code)

    return 