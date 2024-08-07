import numpy as np


def value_decomposition(
    Z,
    X,
    U,
    V,
    T,
    pee,
    pa,
    b,
    theta,
    delta=0.02,
    kappa=2.094215255,
    zeta=1.66e-4 * 1e9,
):
    # Compute change in X
    X_dot = np.diff(X, axis=0)

    # Compute agricultural output
    results_AO = [pa * np.dot(Z[t + 1], theta) / ((1 + delta) ** t) for t in range(T)]
    total_AO = np.sum(results_AO)

    # Compute net transfers
    results_NT = [
        (-b * (kappa * np.sum(Z[t + 1]) - np.sum(X_dot[t])) / ((1 + delta) ** t))
        for t in range(T)
    ]
    total_NT = np.sum(results_NT)

    # Compute forest services
    results_FS = [
        -pee * (kappa * np.sum(Z[t + 1]) - np.sum(X_dot[t])) / ((1 + delta) ** t)
        for t in range(T)
    ]
    total_FS = np.sum(results_FS)

    # Compute adjustment costs
    results_AC = [
        (zeta / 2) * (np.sum(U[t]) + np.sum(V[t])) ** 2 / ((1 + delta) ** t)
        for t in range(T)
    ]
    total_AC = np.sum(results_AC)

    # Compute total net present value
    total_PV = total_AO + total_NT + total_FS - total_AC

    return {
        "pa": pa,
        "pee": pee,
        "b": b,
        "total_AO": total_AO,
        "total_NT": total_NT,
        "total_FS": total_FS,
        "total_AC": total_AC,
        "total_PV": total_PV,
    }


def transfers_decomposition(
    years,
    Z,
    X,
    Z_base,
    X_base,
    b,
    delta=0.02,
    kappa=2.094215255,
):
    # Compute change in X
    X_dot = np.diff(X, axis=0)
    X_dot_base = np.diff(X_base, axis=0)

    # Compute net captured emissions for base case
    results_NCE_base = [-kappa * Z_base[t + 1] + X_dot_base[t] for t in range(years)]
    total_NCE_base = np.sum(results_NCE_base)

    # Compute NCE
    results_NCE = [-kappa * Z[t + 1] + X_dot[t] for t in range(years)]
    total_NCE = np.sum(results_NCE)

    results_NT2 = [
        -b * (kappa * Z[t + 1] - X_dot[t]) / ((1 + delta) ** t) for t in range(years)
    ]
    total_NT2 = np.sum(results_NT2)

    total_EC = total_NT2 / (total_NCE - total_NCE_base)

    return {
        "b": b,
        "net captured emissions": total_NCE,
        "discounted net transfers": total_NT2,
        "discounted effective costs": total_EC,
    }
