import numpy as np
from cmdstanpy import CmdStanModel

from ..sampling import gamma_adj_reg_data, theta_adj_reg_data
from ..services.data_service import load_productivity_reg_data
from ..services.file_service import get_path


def sample(num_sites: int, **stan_kwargs):
    # Load parameter regression data
    (
        site_theta_df,
        site_gamma_df,
        municipal_theta_df,
        municipal_gamma_df,
    ) = load_productivity_reg_data(num_sites)

    # Read model code
    stan_file = get_path("stan_model") / "baseline.stan"
    sampler = CmdStanModel(stan_file=stan_file, force_compile=True)

    # Pack into model data
    model_data = dict(
        S=num_sites,
        pa_2017=44.9736197781184,
        **theta_adj_reg_data(num_sites, site_theta_df),
        **gamma_adj_reg_data(num_sites, site_gamma_df),
        **baseline_hyperparams(municipal_theta_df, "theta"),
        **baseline_hyperparams(municipal_gamma_df, "gamma"),
    )

    # Sampling
    fit = sampler.sample(
        data=model_data,
        fixed_param=True,
        **stan_kwargs,
    )
    return fit


def baseline_hyperparams(municipal_df, var):
    # Drop records with missing data
    municipal_df = municipal_df.dropna()

    if var == "theta":
        # Get theta regression data
        y, X, W = _theta_muni_reg_data(municipal_df)

        # Applying WLS weights
        y = W @ y
        X = W @ X

    elif var == "gamma":
        # Get gamma regression data
        y, X = _gamma_muni_reg_data(municipal_df)
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


def _theta_muni_reg_data(df):
    # Get outcome
    y = df["log_slaughter_value_per_ha"].to_numpy()

    # Normalize weights
    w = df["weights"] / df["weights"].sum()

    # # Get weights pre-multiplication matrix
    W = np.diag(np.sqrt(w))

    # Get regression design matrix and its dimensions
    X = df.iloc[:, :8].to_numpy()

    return y, X, W


def _gamma_muni_reg_data(df):
    # Get outcome
    y = df["log_co2e_ha_2017"].to_numpy()

    # Get regression design matrix and its dimensions
    X = df.iloc[:, :5].to_numpy()
    return y, X
