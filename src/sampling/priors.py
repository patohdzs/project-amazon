import stan
from services.data_service import load_site_data
from services.file_service import stan_model_path
from solvers.stan import _gamma_reg_data, _prior_hyperparams, _theta_reg_data


def sample_priors(model_name: str, num_sites: str):
    # Load sites data
    (
        zbar_2017,
        gamma_vals,
        z_2017,
        forestArea_2017_ha,
        theta_vals,
        gamma_coe_mean,
        theta_coe_mean,
        gamma_coe_vcov,
        theta_coe_vcov,
        theta_data,
        gamma_data,
    ) = load_site_data(num_sites, norm_fac=1e11)

    # Read model code
    with open(stan_model_path(model_name) / "prior.stan") as f:
        model_code = f.read()

    # Get regression data
    y_theta, X_theta, N_theta, K_theta, G_theta = _theta_reg_data(num_sites, theta_data)
    y_gamma, X_gamma, N_gamma, K_gamma, G_gamma = _gamma_reg_data(num_sites, gamma_data)

    # Pack into model data
    model_data = dict(
        S=num_sites,
        K_theta=K_theta,
        K_gamma=K_gamma,
        N_theta=N_theta,
        N_gamma=N_gamma,
        X_theta=X_theta,
        X_gamma=X_gamma,
        G_theta=G_theta,
        G_gamma=G_gamma,
        pa_2017=44.9736197781184,
        **_prior_hyperparams(y_theta, X_theta, "theta"),
        **_prior_hyperparams(y_gamma, X_gamma, "gamma"),
    )

    # Compiling model
    sampler = stan.build(program_code=model_code, data=model_data, random_seed=1)

    # Sampling
    fit = sampler.fixed_param(num_samples=1000)
    return fit
