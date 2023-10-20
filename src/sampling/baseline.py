import numpy as np
import stan
from sampling import gamma_reg_data, theta_reg_data
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
    ) = load_site_data(num_sites)

    # Read model code
    with open(stan_model_path(model_name) / "baseline.stan") as f:
        model_code = f.read()

    # Get regression data
    _, X_theta, N_theta, K_theta, G_theta, _ = theta_reg_data(num_sites, theta_data)
    _, X_gamma, N_gamma, K_gamma, G_gamma = gamma_reg_data(num_sites, gamma_data)

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
        **baseline_hyperparams(num_sites, theta_data, "theta"),
        **baseline_hyperparams(num_sites, gamma_data, "gamma"),
    )

    # Compiling model
    sampler = stan.build(program_code=model_code, data=model_data, random_seed=1)

    # Sampling
    fit = sampler.fixed_param(num_samples=num_samples)
    return fit


def baseline_hyperparams(num_sites, df, var):
    # Drop records with missing data
    df = df.dropna()

    if var == "theta":
        # Get theta regression data
        y, X, _, _, _, W = theta_reg_data(num_sites, df)

        # Applying WLS weights
        y = W @ y
        X = W @ X

    elif var == "gamma":
        # Get gamma regression data
        y, X, _, _, _ = gamma_reg_data(num_sites, df)
    else:
        raise ValueError("Argument `var` should be one of `theta`, `gamma`")

    inv_Q = np.linalg.inv(X.T @ X)
    m = inv_Q @ X.T @ y
    a = (X.shape[0]) / 2
    b = 0.5 * (y.T @ y - m.T @ X.T @ X @ m)
    return {
        f"inv_Q_{var}": inv_Q,
        f"m_{var}": m,
        f"a_{var}": a,
        f"b_{var}": b,
    }
