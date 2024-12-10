import matplotlib.pyplot as plt
from pysrc.services.data_service import load_productivity_params
from pysrc.optimization import solve_planner_problem
from pysrc.services.data_service import load_site_data
from pysrc.analysis import value_decomposition

# Model hyperparameters
solver = "gurobi"
pee = 6.6
pa = 41.11
num_sites = 78
T = 200
b=0

# Load site data
(zbar_2017, z_2017, forest_area_2017) = load_site_data(num_sites)



(theta, gamma) = load_productivity_params(num_sites)
# Computing carbon absorbed in start period
x_2017 = gamma * forest_area_2017



results= solve_planner_problem(
        x0=x_2017,
        z0=z_2017,
        zbar=zbar_2017,
        gamma=gamma,
        theta=theta,
        time_horizon=T,
        price_cattle=pa,
        price_emissions=pee + b,
    )
    

print(
    "result",
    value_decomposition(
        solution=results,
        T=T,
        pee=pee,
        pa=pa,
        b=b,
        theta=theta,
    ),
)