import pandas as pd

from pysrc.analysis import value_decomposition
from pysrc.optimization.gurobi import solve_planner_problem
from pysrc.sampling import baseline
from pysrc.services.data_service import load_site_data

# ## Model scenario


solver = "gurobi"
pee = 7.1
pa = 41.11
sitenum = 78
T = 200
b = 25

# Load site data
(
    zbar_2017,
    z_2017,
    forest_area_2017,
    _,
    _,
    _,
    _,
) = load_site_data(sitenum)

# Set productivity parameters using baseline mean
baseline_fit = baseline.sample(
    num_sites=sitenum,
    iter_sampling=10**4,
    chains=5,
    seed=1,
)

theta = baseline_fit.stan_variable("theta").mean(axis=0)
gamma = baseline_fit.stan_variable("gamma").mean(axis=0)


base_parameter = pd.DataFrame({"theta": theta.round(2), "gamma": gamma.round(2)})

# Save the DataFrame to a CSV file
base_parameter.to_csv("theta_gamma.csv", index=False)

# Computing carbon absorbed in start period
x0_vals = gamma * forest_area_2017


results = solve_planner_problem(
    T=T,
    theta=theta,
    gamma=gamma,
    x0=x0_vals,
    zbar=zbar_2017,
    z0=z_2017,
    pe=pee + b,
    pa=pa,
)


print(
    "result",
    value_decomposition(
        Z=results["Z"],
        X=results["X"],
        U=results["U"],
        V=results["V"],
        T=T,
        pee=pee,
        pa=pa,
        b=b,
        theta=theta,
    ),
)
