from pysrc.analysis import value_decomposition
from pysrc.optimization.gurobi import solve_planner_problem
from pysrc.sampling import baseline
from pysrc.services.data_service import load_site_data

# Model hyperparameters
solver = "gurobi"
pee = 7.6
pa = 41.11
sitenum = 1043
T = 200
b = 25

# Load site data
(zbar_2017, z_2017, forest_area_2017) = load_site_data(sitenum)

# Set productivity parameters using baseline mean
baseline_fit = baseline.sample(
    num_sites=sitenum,
    iter_sampling=10**4,
    chains=5,
    seed=1,
)

theta = baseline_fit.stan_variable("theta").mean(axis=0)
gamma = baseline_fit.stan_variable("gamma").mean(axis=0)


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
    )
)
