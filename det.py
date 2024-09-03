from pysrc.analysis import value_decomposition
from pysrc.optimization.gams import solve_planner_problem
from pysrc.sampling import baseline
from pysrc.services.data_service import load_site_data

# Model hyperparameters
solver = "gams"
pee = 7.6
pa = 41.11
sitenum = 1043
T = 200
b = 10

# Load site data
(zbar_2017, z_2017, forest_area_2017,_,_,_,_) = load_site_data(sitenum)

# Set productivity parameters using baseline mean
baseline_fit = baseline.sample(
    num_sites=sitenum,
    iter_sampling=10**3,
    chains=5,
    seed=1,
)

theta = baseline_fit.stan_variable("theta").mean(axis=0)
gamma = baseline_fit.stan_variable("gamma").mean(axis=0)

# Computing carbon absorbed in start period
x_2017 = gamma * forest_area_2017

# Solve planner problem
results = solve_planner_problem(
    x0=x_2017,
    z0=z_2017,
    zbar=zbar_2017,
    gamma=gamma,
    theta=theta,
    T=T,
    pe=pee + b,
    pa=pa,
)



print(
    value_decomposition(
        Z=results['Z'],
        X=results['X'],
        U=results['U'],
        V=results['V'],
        T=T,
        pee=pee,
        pa=pa,
        b=b,
        theta=theta,
    )
)