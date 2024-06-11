
import geopandas as gpd
import os
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
import numpy as np
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
import pandas as pd
from pysrc.services.file_service import get_path

def spatial_allocation(opt='gams',
                       pa=41.11,
                       pe_hmc=7.1,
                       pe_det=5.3,
                       b=0):
    pe_hmc+=b
    pe_det+=b
    num_sites=78
    # Load the data
    clean_amazonBiome = gpd.read_file(str(get_path("data"))+'/hmc/map.geojson')
    id=gpd.read_file(str(get_path("data"))+"/hmc/id_78.geojson").to_crs(epsg=4326)
    id_hmc = id.copy()
    id_det = id.copy()
    id_hmc_half=id.copy()
    id_det_half=id.copy()


    result_folder = os.path.join(str(get_path("output")), "optimization")
    output_folder = str(get_path("output"))+"/figures/maps/"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    z_hmc= pd.read_csv(result_folder+'/hmc/'+opt+f"/78sites/pa_{pa}/"+f"pe_{pe_hmc}/amazon_data_z.dat", delimiter='\t')
    z_hmc=z_hmc.drop('T/R ', axis=1).to_numpy().T*1e9
    z_det= pd.read_csv(result_folder+'/det/'+opt+f"/78sites/pa_{pa}/"+f"pe_{pe_det}/amazon_data_z.dat", delimiter='\t')
    z_det=z_det.drop('T/R ', axis=1).to_numpy().T*1e9


    u_hmc= pd.read_csv(result_folder+'/hmc/'+opt+f"/78sites/pa_{pa}/"+f"pe_{pe_hmc}/amazon_data_u.dat", delimiter='\t')
    u_hmc=u_hmc.drop('T/R ', axis=1).to_numpy().T[:,:50]*1e9
    v_hmc= pd.read_csv(result_folder+'/hmc/'+opt+f"/78sites/pa_{pa}/"+f"pe_{pe_hmc}/amazon_data_v.dat", delimiter='\t')
    v_hmc=v_hmc.drop('T/R ', axis=1).to_numpy().T[:,:50]*1e9

    u_det= pd.read_csv(result_folder+'/det/'+opt+f"/78sites/pa_{pa}/"+f"pe_{pe_det}/amazon_data_u.dat", delimiter='\t')
    u_det=u_det.drop('T/R ', axis=1).to_numpy().T[:,:50]*1e9
    v_det= pd.read_csv(result_folder+'/det/'+opt+f"/78sites/pa_{pa}/"+f"pe_{pe_det}/amazon_data_v.dat", delimiter='\t')
    v_det=v_det.drop('T/R ', axis=1).to_numpy().T[:,:50]*1e9
    hmc_max = np.maximum(np.abs(u_hmc), np.abs(v_hmc))
    det_max = np.maximum(np.abs(u_det), np.abs(v_det))
    
    indicator_hmc = np.full(z_hmc.shape[0], 2)
    indicator_hmc[z_hmc[:, 49] > z_hmc[:, 0]] = 0
    indicator_hmc[z_hmc[:, 49] < z_hmc[:, 0]] = 1
    indicator_hmc[np.abs(z_hmc[:,49]-z_hmc[:,0])<0.05*z_hmc[:,0]]=2
    
    indicator_det = np.full(z_det.shape[0], 2)
    indicator_det[z_det[:, 49] > z_det[:, 0]] = 0
    indicator_det[z_det[:, 49] < z_det[:, 0]] = 1
    indicator_det[np.abs(z_det[:,49]-z_det[:,0])<0.05*z_det[:,0]]=2
    
    positions_hmc = []

    for i in range(hmc_max.shape[0]):
        series = hmc_max[i, :50]
        found = False  
        for j in range(0, len(series)):
            if series[j] ==np.max(series) and series[j]>1e-4:
                positions_hmc.append(j+1)
                found = True
                break  
        
        if not found:  
            positions_hmc.append(111)


    positions_det = []

    for i in range(det_max.shape[0]):
        series = det_max[i, :50]
        found = False 
        for j in range(0, len(series)):
            if series[j] ==np.max(series) and series[j]>1e-4:

                positions_det.append((j+1))
                found = True
                break  
        if not found:  
            positions_det.append(111)



    id_hmc['id']=positions_hmc
    id_det['id']=positions_det
    id_index=id.copy()

    positions_index=[]
    for i in range(1,79):
        positions_index.append(i)
    id_index['id']=positions_index



    fig, ax = plt.subplots()

    id_hmc.boundary.plot(edgecolor='black', facecolor='none', linewidth=0.5, ax=ax)
    clean_amazonBiome.boundary.plot(color="black", linewidth=1.2, ax=ax)

    shaded_areas = [
        (-53.366934735915706, -12.921187383718369),
        (-60.64328853728383, -12.921187383718369)
    ]

    for x, y, hmc_id, det_id,indicator in zip(id_hmc.geometry.centroid.x, id_hmc.geometry.centroid.y, id_hmc['id'], id_det['id'],indicator_hmc):

        if indicator ==0:
            color = 'red'
        elif indicator==1:
            color = 'green'
        elif indicator==2:
            hmc_id=""  
            
        ax.text(x, y, str(hmc_id), fontsize=14, ha='center', va='center', fontweight='bold', color=color)

    for x, y in shaded_areas:
        rect = Rectangle((x - 1.2, y - 1.2), 2.4, 2.4, color='grey', alpha=0.3, linewidth=0, fill=True)
        ax.add_patch(rect)

    ax.axis('off') 

    plt.tight_layout()
    fig.savefig(output_folder+f'map_78site_ambiguity_b{b}_pehmc_{pe_hmc}_pedet_{pe_det}.png', dpi=300)
    plt.show()



    fig, ax = plt.subplots()

    id_det.boundary.plot(edgecolor='black', facecolor='none', linewidth=0.5, ax=ax)

    clean_amazonBiome.boundary.plot(color="black", linewidth=1.2, ax=ax)

    for x, y, hmc_id, det_id,indicator in zip(id_det.geometry.centroid.x, id_det.geometry.centroid.y, id_hmc['id'], id_det['id'],indicator_det):
        if indicator ==0:
            color = 'red'
        elif indicator==1:
            color = 'green'
        elif indicator==2:
            det_id=""  

        ax.text(x, y, str(det_id), fontsize=14, ha='center', va='center', fontweight='bold', color=color)

    for x, y in shaded_areas:
        rect = Rectangle((x - 1.2, y - 1.2), 2.4, 2.4, color='grey', alpha=0.3, linewidth=0, fill=True)
        ax.add_patch(rect)

    ax.axis('off') 

    plt.tight_layout()
    fig.savefig(output_folder+f'map_78site_neutral_b{b}_pehmc_{pe_hmc}_pedet_{pe_det}.png', dpi=300)
    plt.show()

    return

