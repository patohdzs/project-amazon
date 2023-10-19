import stan
from sampling.stan import _gamma_reg_data, _prior_hyperparams, _theta_reg_data
from services.data_service import load_site_data
from services.file_service import stan_model_path


def sample(model_name: str, num_samples: int, num_sites: int):
    # Load sites data
    (
        _,
        _,
        _,
        _,
        _,
        theta_data,
        gamma_data,
    ) = load_site_data(num_sites, norm_fac=1e11)

    # Read model code
    with open(stan_model_path(model_name) / "prior.stan") as f:
        model_code = f.read()

    # Get regression data
    _, X_theta, N_theta, K_theta, G_theta, _ = _theta_reg_data(num_sites, theta_data)
    _, X_gamma, N_gamma, K_gamma, G_gamma = _gamma_reg_data(num_sites, gamma_data)

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
        **_prior_hyperparams(num_sites, theta_data, "theta"),
        **_prior_hyperparams(num_sites, gamma_data, "gamma"),
    )

    # Compiling model
    sampler = stan.build(program_code=model_code, data=model_data, random_seed=1)

    # Sampling
    fit = sampler.fixed_param(num_samples=num_samples)
    return fit
