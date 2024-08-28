from pysrc.sampling import baseline
from pysrc.services.data_service import load_site_data
import pandas as pd
num_sites=1043

baseline_fit = baseline.sample(
    num_sites=num_sites, iter_sampling=10**4, chains=5, seed=1
)
# theta = baseline_fit.stan_variable("beta_theta")
gamma = baseline_fit.stan_variable("beta_gamma")
# gamma = baseline_fit.stan_variable("gamma").mean(axis=0)
gamma_df = pd.DataFrame(gamma)

# gamma_df = pd.DataFrame(gamma,columns=['intercept','log_hist precip','log_hist_temp','lat','lon','latlon'])

# Save to CSV files
gamma_df.to_csv('gamma.csv', index=False)


# print("theta",theta)
# print("gamma",gamma)


# (
#     _,
#     _,
#     _,
#     site_theta_df,
#     site_gamma_df,
#     municipal_theta_df,
#     municipal_gamma_df,
# ) = load_site_data(78)

# print(site_gamma_df)