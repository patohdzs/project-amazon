import numpy as np


def value_decomposition(
    Z,
    X,
    U,
    V,
    T,
    theta,
    delta=0.02,
    kappa=2.094215255,
    zeta=1.66e-4 * 1e9,
    pee=20.76,
    b=0,
    pa=44.75,
):
    # Compute change in X
    X_dot = np.diff(X, axis=0)

    # Compute agricultural output
    results_AO = [
        pa * np.dot(Z[t], theta[0]) / ((1 + delta) ** (t + 1)) for t in range(T)
    ]
    total_AO = np.sum(results_AO)

    # Compute net transfers
    results_NT = [
        (-b * (kappa * np.sum(Z[t]) - X_dot[t]) / ((1 + delta) ** (t + 1)))
        for t in range(T)
    ]
    total_NT = np.sum(results_NT)

    # Compute climate services
    results_CS = [
        -pee * (kappa * np.sum(Z[t]) - X_dot[t]) / ((1 + delta) ** (t + 1))
        for t in range(T)
    ]
    total_CS = np.sum(results_CS)

    # Compute adjustment costs
    results_AC = [
        (zeta / 2) * (np.sum(U[t]) + np.sum(V[t])) ** 2 / ((1 + delta) ** (t + 1))
        for t in range(T)
    ]
    total_AC = np.sum(results_AC)

    # Compute total net present value
    total_PV = total_AO + total_NT + total_CS - total_AC

    return {
        "pa": pa,
        "pee": pee,
        "b": b,
        "total_AO": total_AO,
        "total_NT": total_NT,
        "total_CS": total_CS,
        "total_AC": total_AC,
        "total_PV": total_PV,
    }
