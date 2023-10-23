import numpy as np
import stan
from sampling import gamma_adj_reg_data, theta_adj_reg_data
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
        site_theta_df,
        site_gamma_df,
        municipal_theta_df,
        municipal_gamma_df,
    ) = load_site_data(num_sites)

    # Read model code
    with open(stan_model_path(model_name) / "baseline.stan") as f:
        model_code = f.read()

    # Get regression data
    _, X_theta, N_theta, K_theta, G_theta, _ = theta_adj_reg_data(
        num_sites, site_theta_df
    )
    _, X_gamma, N_gamma, K_gamma, G_gamma = gamma_adj_reg_data(num_sites, site_gamma_df)

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
        **baseline_hyperparams(municipal_theta_df, "theta"),
        **baseline_hyperparams(municipal_gamma_df, "gamma"),
    )

    # Compiling model
    sampler = stan.build(program_code=model_code, data=model_data, random_seed=1)

    # Sampling
    fit = sampler.fixed_param(num_samples=num_samples)
    return fit


def baseline_hyperparams(municipal_df, var):
    # Drop records with missing data
    municipal_df = municipal_df.dropna()

    if var == "theta":
        # Get theta regression data
        y, X, W = theta_baseline_reg_data(municipal_df)

        # Applying WLS weights
        y = W @ y
        X = W @ X

    elif var == "gamma":
        # Get gamma regression data
        y, X = gamma_baseline_reg_data(municipal_df)
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


def theta_baseline_reg_data(municipal_theta_df):
    # Get outcome
    y = municipal_theta_df["log_cattleSlaughter_valuePerHa_2017"].to_numpy()

    # Get weights matrix
    W = np.diag(np.sqrt(municipal_theta_df["weights"]))

    # Get regression design matrix and its dimensions
    X = municipal_theta_df.iloc[:, :8].to_numpy()

    return y, X, W


def gamma_baseline_reg_data(municipal_gamma_df):
    # Get outcome
    y = municipal_gamma_df["log_co2e_ha_2017"].to_numpy()

    # Get regression design matrix and its dimensions
    X = municipal_gamma_df.iloc[:, :5].to_numpy()
    return y, X
