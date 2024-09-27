import numpy as np
import pandas as pd
import geopandas as gpd
from scipy.stats import invgamma, norm,gamma
import matplotlib.pyplot as plt
from pysrc.services.file_service import get_path

data_dir = get_path("data", "calibration", "hmc")
municipal_gamma_df = gpd.read_file(data_dir/"gamma_reg_site_1043.geojson")

X = municipal_gamma_df.iloc[:, 0:6].values 
Y = municipal_gamma_df.iloc[:, 6].values   
site_ids = municipal_gamma_df.iloc[:, 7].values -1   

np.random.seed(45)



def gibbs_sampling_with_data(X, Y, site_ids, n_iterations=500000, beta_prior_mean=None, beta_prior_cov=None, 
                             eta_prior_shape=2, eta_prior_scale=1, nu_prior_shape=2, nu_prior_scale=1):
    n_data, n_features = X.shape
    n_sites = len(np.unique(site_ids))
    
    # Initialize priors
    if beta_prior_mean is None:
        beta_prior_mean = np.zeros(n_features)
    if beta_prior_cov is None:
        beta_prior_cov = np.eye(n_features)   # Large prior variance
    
    # Initialize Gibbs sampling arrays
    beta_samples = np.zeros((n_iterations, n_features))
    eta_samples = np.zeros(n_iterations)
    nu_samples = np.zeros(n_iterations)
    V_samples = np.zeros((n_iterations, n_sites))
    

    
    
    # Initial values based on prior
    beta_current = np.random.multivariate_normal(beta_prior_mean, beta_prior_cov)
    eta_current = gamma.rvs(a=eta_prior_shape, scale=1/eta_prior_scale)
    nu_current = gamma.rvs(a=nu_prior_shape, scale=1/nu_prior_scale)
    V_current = np.zeros(n_sites)+0.01
    
    # Gibbs sampler
    for i in range(n_iterations):
        # Step i: Draw beta conditioned on V, eta, and data
        precision_beta = eta_current * (X.T @ X)
        mean_beta = np.linalg.inv(precision_beta) @ ( eta_current * X.T @ (Y - V_current[site_ids]))
        beta_current = np.random.multivariate_normal(mean_beta, np.linalg.inv(precision_beta))
        beta_samples[i, :] = beta_current
        
        # Step ii: Draw eta conditioned on V and data
        residuals = Y - X @ beta_current - V_current[site_ids]
        eta_shape_post =    n_data / 2
        eta_scale_post =  0.5 * np.sum(residuals**2)
        eta_current = gamma.rvs(a=eta_shape_post, scale=1/eta_scale_post)
        eta_samples[i] = eta_current
        
        # Step iii: Draw nu conditioned on V
        nu_shape_post =  n_sites / 2
        nu_scale_post =  0.5 * np.sum(V_current**2)
        nu_current = gamma.rvs(a=nu_shape_post, scale=1/nu_scale_post)
        nu_samples[i] = nu_current
        

        # Step iv: Draw V_j conditioned on beta, nu, eta, and data
        for j in range(n_sites):
            site_mask = site_ids == j
            n_j = np.sum(site_mask)
            V_bar_j = np.mean(Y[site_mask] - X[site_mask] @ beta_current)
            mean_V_j = (n_j * eta_current * V_bar_j) / (n_j * eta_current + nu_current)
            precision_V_j = n_j * eta_current + nu_current
            V_current[j] = np.random.normal(mean_V_j, np.sqrt(1 / precision_V_j))
        V_samples[i, :] = V_current
    
    
    beta_samples_final = beta_samples[200000:,:]
    eta_samples_final = eta_samples[200000:]
    nu_samples_final = nu_samples[200000:]
    V_samples_final = V_samples[200000:,:]
    
    # Posterior means
    beta_posterior_mean = beta_samples_final.mean(axis=0)
    eta_posterior_mean = eta_samples_final.mean()
    nu_posterior_mean = nu_samples_final.mean()
    
    return beta_posterior_mean, eta_posterior_mean, nu_posterior_mean, beta_samples_final, eta_samples_final, nu_samples_final, V_samples_final


beta_posterior_mean, eta_posterior_mean, nu_posterior_mean, beta_samples, eta_samples, nu_samples, V_samples = gibbs_sampling_with_data(X, Y, site_ids)


print(f"Posterior mean of beta: {beta_posterior_mean}")
print(f"Posterior mean of eta: {eta_posterior_mean}")
print(f"Posterior mean of nu: {nu_posterior_mean}")



beta_posterior_mean = np.mean(beta_samples,axis=0)
beta_posterior_variance = np.cov(beta_samples, rowvar=False, ddof=0)  

eta_posterior_mean = np.mean(eta_samples)
eta_posterior_variance = np.var(eta_samples, ddof=0)  

nu_posterior_mean = np.mean(nu_samples)
nu_posterior_variance = np.var(nu_samples, ddof=0)  

V_posterior_mean = np.mean(V_samples,axis=0)
V_posterior_variance = np.var(V_samples, ddof=0,axis=0)



data = []


for i in range(len(beta_posterior_mean)):
    data.append({'Parameter': f'beta_{i+1}', 'Posterior Mean': beta_posterior_mean[i], 'Posterior Variance': None})
    

for i in range(len(beta_posterior_mean)):
    for j in range(len(beta_posterior_mean)):
        data.append({'Parameter': f'var_beta_{i+1}_{j+1}', 'Posterior Mean': None, 'Posterior Variance': beta_posterior_variance[i, j]})


data.append({'Parameter': 'eta', 'Posterior Mean': eta_posterior_mean, 'Posterior Variance': eta_posterior_variance})


data.append({'Parameter': 'nu', 'Posterior Mean': nu_posterior_mean, 'Posterior Variance': nu_posterior_variance})


for i in range(len(V_posterior_mean)):
    data.append({'Parameter': f'V_{i+1}', 'Posterior Mean': V_posterior_mean[i], 'Posterior Variance': V_posterior_variance[i]})


df = pd.DataFrame(data)


df.to_csv(data_dir/'gamma_posterior_means_and_variances.csv', index=False)

gamma_fit_df_1043 = gpd.read_file(data_dir/"gamma_data_site_1043.geojson")
X_fit = gamma_fit_df_1043.iloc[:, 0:6].values  
id_fit = gamma_fit_df_1043.iloc[:, 7].values 
gamma_fit_df_1043 = np.mean(np.exp(X_fit @ beta_samples.T + V_samples[:, id_fit-1].T),axis=1)
gamma_fit_df_1043_df = pd.DataFrame(gamma_fit_df_1043, columns=['gamma_fit'])
gamma_fit_df_1043_df.to_csv(data_dir/'gamma_fit_1043.csv', index=False)

gamma_fit_df_78 = gpd.read_file(data_dir/"gamma_data_site_78.geojson")
X_fit = gamma_fit_df_78.iloc[:, 0:6].values  # Columns 1 to 6 as X
id_fit = gamma_fit_df_78.iloc[:, 7].values  # Columns 1 to 6 as X
gamma_fit_df_78 = np.mean(np.exp(X_fit @ beta_samples.T + V_samples[:,].T),axis=1)
gamma_fit_df_78_df = pd.DataFrame(gamma_fit_df_78, columns=['gamma_fit'])
gamma_fit_df_78_df.to_csv(data_dir/'gamma_fit_78.csv', index=False)


print("samping done")