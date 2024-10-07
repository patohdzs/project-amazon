from pysrc.analysis import value_decomposition
from pysrc.optimization import solve_planner_problem,vectorize_trajectories
from pysrc.services.data_service import load_site_data,load_reg_data,load_productivity_params

# Model hyperparameters
solver = "gams"
pee = 7.6
pa = 41.11
sitenum = 78
T = 200
b = 10

# Load site data
(zbar_2017, z_2017, forest_area_2017) = load_site_data(sitenum)


(theta, gamma) = load_productivity_params(78)
# Computing carbon absorbed in start period
x_2017 = gamma * forest_area_2017

# Solve planner problem
results = solve_planner_problem(
    x0=x_2017,
    z0=z_2017,
    zbar=zbar_2017,
    gamma=gamma,
    theta=theta,
    price_emissions=pee + b,
    price_cattle=pa,
)



print(
    value_decomposition(
        Z=results.Z,
        X=results.X,
        U=results.U,
        V=results.V,
        T=T,
        pee=pee,
        pa=pa,
        b=b,
        theta=theta,
    )
)