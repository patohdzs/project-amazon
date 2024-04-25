import os

import numpy as np
import pandas as pd

from pysrc.analysis import value_decomposition

pa = 41.11
pee = 7.1
b = [0, 10, 15, 20, 25]
pe = [pee + bi for bi in b]
kappa = 2.094215255
zeta = 1.66e-4 * 1e11


output_directory = os.getcwd() + "\\"
baseline_directory = os.getcwd() + f"\\gams\\p_a_{pa}_p_e_{pe[0]}\\"
gams_directory = os.getcwd() + "\\gams\\"


def read_file(result_directory):
    dfz = pd.read_csv(result_directory + "\\amazon_data_z.dat", delimiter="\t")
    dfz = dfz.drop("T/R ", axis=1)
    dfz_np = dfz.to_numpy()

    dft = pd.read_csv(result_directory + "\\theta.txt", header=None)
    dft_np = dft.to_numpy().T

    dfx = pd.read_csv(result_directory + "\\amazon_data_x.dat", delimiter="\t")
    dfx = dfx.drop("T   ", axis=1)
    dfx_np = dfx.to_numpy()
    dfxdot = np.diff(dfx_np, axis=0)

    dfu = pd.read_csv(result_directory + "\\amazon_data_u.dat", delimiter="\t")
    dfu = dfu.drop("T/R ", axis=1)
    dfu_np = dfu.to_numpy()

    dfv = pd.read_csv(result_directory + "\\amazon_data_v.dat", delimiter="\t")
    dfv = dfv.drop("T/R ", axis=1)
    dfv_np = dfv.to_numpy()

    return (dfz_np, dft_np, dfxdot, dfu_np, dfv_np)


def table_one(
    Zs,
    Xs,
    Us,
    Vs,
    bs,
    T,
    pee,
    pa,
    theta,
    delta=0.02,
    kappa=2.094215255,
    zeta=1.66e-4 * 1e9,
):
    results = [
        value_decomposition(
            Z=Z,
            X=X,
            U=U,
            V=V,
            T=T,
            pee=pee,
            pa=pa,
            b=b,
            theta=theta,
            delta=delta,
            kappa=kappa,
            zeta=zeta,
        )
        for Z, X, U, V, b in zip(Zs, Xs, Us, Vs, bs)
    ]

    table_df = pd.DataFrame(results)
    return table_df.to_latex(index=False)


def table_two():
    dfz = pd.read_csv(baseline_directory + "amazon_data_z.dat", delimiter="\t")
    dfz = dfz.drop("T/R ", axis=1)
    dfz_np = dfz.to_numpy()

    dfx = pd.read_csv(baseline_directory + "amazon_data_x.dat", delimiter="\t")
    dfx = dfx.drop("T   ", axis=1)
    dfx_np = dfx.to_numpy()
    dfxdot = np.diff(dfx_np, axis=0)

    results_NCE_base = []
    for i in range(30):
        result_NCE_base = -kappa * np.sum(dfz_np[i + 1]) + dfxdot[i]
        results_NCE_base.append(result_NCE_base)
    total_NCE_base = np.sum(results_NCE_base) * 100

    results_table2 = []
    for j in range(5):
        order = j

        result_directory = gams_directory + f"p_a_{pa}_p_e_{pe[order]}"

        (dfz_np, dft_np, dfxdot, dfu_np, dfv_np) = read_file(result_directory)

        results_NCE = []
        for i in range(30):
            result_NCE = -kappa * np.sum(dfz_np[i + 1]) + dfxdot[i]
            results_NCE.append(result_NCE)
        total_NCE = np.sum(results_NCE) * 100

        results_NT2 = []
        for i in range(30):
            result_NT2 = (
                -b[order]
                * (kappa * np.sum(dfz_np[i + 1]) - dfxdot[i])
                / ((1 + 0.02) ** (i))
            )

            results_NT2.append(result_NT2)
        total_NT2 = np.sum(results_NT2)

        total_EC = total_NT2 / (total_NCE - total_NCE_base) * 100

        total_NCE = total_NCE.round(2)
        total_NT2 = total_NT2.round(2)
        total_EC = total_EC.round(2)
        iteration_results = {
            "pe": pe[order],
            "b": b[order],
            "total_NCE": total_NCE,
            "total_NT2": total_NT2,
            "total_EC": total_EC,
        }

        results_table2.append(iteration_results)

    results_df2 = pd.DataFrame(results_table2)
    return results_df2.to_latex(index=False)
