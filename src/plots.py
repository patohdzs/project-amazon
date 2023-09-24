import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def traceplot_abs_error(results: dict, plots_dir: Path) -> None:
    # Make plot fig and axis
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    # Lineplot
    plt.plot(results["abs_error_tracker"], label=r"Absolute Error")

    # Labels, title and legend
    plt.xlabel("Iteration")
    plt.ylabel(r"Absolute Error")
    plt.title(r"Trace Plot of Absolute Error")
    legend = plt.legend(bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0)
    fig.tight_layout()
    plt.subplots_adjust(right=0.7)

    # Save figure
    fig.savefig(
        plots_dir / "abs_error.png",
        bbox_extra_artists=(legend,),
        bbox_inches="tight",
        dpi=100,
    )
    plt.close()


def traceplot_pct_error(results: dict, plots_dir: Path) -> None:
    # Make plot fig and axis
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    plt.plot(results["percentage_error_tracker"], label=r"Proportional Error")
    plt.xlabel("Iteration")
    plt.ylabel(r"Proportional Error")
    plt.title(r"Trace Plot of Proportional Error")
    legend = plt.legend(bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0)
    fig.tight_layout()
    plt.subplots_adjust(right=0.7)
    fig.savefig(
        plots_dir / "pro_error.png",
        bbox_extra_artists=(legend,),
        bbox_inches="tight",
        dpi=100,
    )
    plt.close()


def traceplot_params_pct_error(results: dict, plots_dir: Path, K=8) -> None:
    # Find pct change
    pct_change = (
        abs(np.diff(results["uncertain_vals_tracker"], axis=0))
        / results["uncertain_vals_tracker"][:-1, :]
    )

    # Create mpl base
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    for i in range(K):
        axes[0].plot(pct_change[:, i], label=f"Theta_coef {i+1}")

    axes[0].set_title("Theta_coef Proportional Error")
    for i in range(K, pct_change.shape[1]):
        axes[1].plot(pct_change[:, i], label=f"Gamma_coef {i+1-K}")

    axes[1].set_title("Gamma_coef Proportional Error")
    y_min = np.min(pct_change)
    y_max = np.max(pct_change)
    axes[0].set_ylim([y_min, y_max])
    axes[1].set_ylim([y_min, y_max])
    fig.tight_layout()

    fig.savefig(
        plots_dir / "site_pro_error.png",
        bbox_extra_artists=(plt.legend(),),
        bbox_inches="tight",
        dpi=100,
    )
    plt.close()


def traceplot_gamma_coefs(results: dict, plots_dir: Path, K=8) -> None:
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    for j in range(K):
        plt.plot(
            results["uncertain_vals_tracker"][:, j],
            label=r"$\beta_\gamma^{%d}$" % (j + 1),
        )
    plt.xlabel("Iteration")
    plt.ylabel(r"$\beta_\gamma$")
    plt.title(r"Trace Plot of $\beta_\gamma$")
    legend = plt.legend(bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0)
    fig.tight_layout()
    plt.subplots_adjust(right=0.7)
    fig.savefig(
        plots_dir / "gamma_coef.png",
        bbox_extra_artists=(legend,),
        bbox_inches="tight",
        dpi=100,
    )
    plt.close()


def traceplot_theta_coefs(results: dict, plots_dir: Path, K=8) -> None:
    fig, axes = plt.subplots(1, 1, figsize=(8, 6))
    for j in range(K, 13):
        plt.plot(
            results["uncertain_vals_tracker"][:, j],
            label=r"$\beta_\theta^{%d}$" % (j + 1),
        )
    plt.xlabel("Iteration")
    plt.ylabel(r"$\beta_\theta$")
    plt.title(r"Trace Plot of $\beta_\theta$")
    legend = plt.legend(bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0)
    fig.tight_layout()
    plt.subplots_adjust(right=0.7)
    fig.savefig(
        plots_dir / "theta_coef.png",
        bbox_extra_artists=(legend,),
        bbox_inches="tight",
        dpi=100,
    )
    plt.close()


def Z_trajectory(results: dict, plots_dir: Path) -> None:
    for site in range(results["size"]):
        path = plots_dir / "Z_trajectory" / f"site_{site}"
        if not os.path.exists(path):
            os.makedirs(path)

        for i, Z in enumerate(results["sol_val_X_tracker"]):
            fig, axes = plt.subplots(1, 1, figsize=(8, 6))
            plt.plot(
                Z[site, :],
                label=r"$Z_{site_%d, iter_%d}$" % (site, i),
            )
            plt.xlabel("Iteration")
            plt.ylabel(r"$Z$")
            plt.title(f"Trajectory of $Z_{site}$")
            legend = plt.legend(
                bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0
            )
            fig.tight_layout()
            plt.subplots_adjust(right=0.7)
            fig.savefig(
                path / f"Z_site_{site}_iter_{i}.png",
                bbox_extra_artists=(legend,),
                bbox_inches="tight",
                dpi=100,
            )
            plt.close()


def X_trajectory(results: dict, plots_dir: Path) -> None:
    path = plots_dir / "X_trajectory"
    if not os.path.exists(path):
        os.makedirs(path)

    for i, X in enumerate(results["sol_val_X_tracker"]):
        fig, axes = plt.subplots(1, 1, figsize=(8, 6))
        plt.plot(X[-2, :], label=r"$X$")
        plt.xlabel("Time")
        plt.ylabel(r"$X$")
        plt.title(r"Trajectory of $X$")
        legend = plt.legend(
            bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0
        )
        fig.tight_layout()
        plt.subplots_adjust(right=0.7)
        fig.savefig(
            path / f"X_iter_{i}.png",
            bbox_extra_artists=(legend,),
            bbox_inches="tight",
            dpi=100,
        )
        plt.close()


def delta_Z_trajectory(results: dict, plots_dir: Path) -> None:
    for site in range(results["size"]):
        path = plots_dir / "delta_Z_trajectory" / f"site_{site}"
        if not os.path.exists(path):
            os.makedirs(path)

        for i, Z in enumerate(results["sol_val_Z_tracker"]):
            fig, axes = plt.subplots(1, 1, figsize=(8, 6))
            plt.plot(
                Z[site, :],
                label=r"$\Delta Z_{site_%d, iter_%d}$" % (site, i),
            )
            plt.xlabel("Time")
            plt.ylabel(r"$U - V$")
            plt.title(r"Trajectory of $\Delta Z$")
            legend = plt.legend(
                bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0
            )
            fig.tight_layout()
            plt.subplots_adjust(right=0.7)
            fig.savefig(
                path / f"delta_Z_site_{site}_iter_{i}.png",
                bbox_extra_artists=(legend,),
                bbox_inches="tight",
                dpi=100,
            )
            plt.close()


def adj_costs_trajectory(results: dict, plots_dir: Path) -> None:
    path = plots_dir / "Ua_trajectory"
    if not os.path.exists(path):
        os.makedirs(path)

    for i, Ua in enumerate(results["sol_val_Ua_tracker"]):
        fig, axes = plt.subplots(1, 1, figsize=(8, 6))
        plt.plot(Ua, label=r"$Ua$")

        plt.xlabel("Time")
        plt.ylabel(r"$Ua$")
        plt.title(r"Trajectory of Ua")
        legend = plt.legend(
            bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0
        )
        fig.tight_layout()
        plt.subplots_adjust(right=0.7)
        fig.savefig(
            path / f"Ua_iter_{i}.png",
            bbox_extra_artists=(legend,),
            bbox_inches="tight",
            dpi=100,
        )
        plt.close()


def V_trajectory(results: dict, plots_dir: Path) -> None:
    for site in range(results["size"]):
        path = plots_dir / "V_trajectory" / f"site_{site}"
        if not os.path.exists(path):
            os.makedirs(path)

        for i, V in enumerate(results["sol_val_Um_tracker"]):
            fig, axes = plt.subplots(1, 1, figsize=(8, 6))
            plt.plot(V[site, :], label=r"site_%d_iter_%d" % (site, i))
            plt.xlabel("Time")
            plt.ylabel(r"$V$")
            plt.title(r"Trajectory of V")
            legend = plt.legend(
                bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0
            )
            fig.tight_layout()
            plt.subplots_adjust(right=0.7)
            fig.savefig(
                path / f"V_site_{site}_iter_{i}.png",
                bbox_extra_artists=(legend,),
                bbox_inches="tight",
                dpi=100,
            )
            plt.close()


def U_trajectory(results: dict, plots_dir: Path) -> None:
    for site in range(results["size"]):
        path = plots_dir / "U_trajectory" / f"site_{site}"
        if not os.path.exists(path):
            os.makedirs(path)

        for i, U in enumerate(results["sol_val_Up_tracker"]):
            fig, axes = plt.subplots(1, 1, figsize=(8, 6))
            plt.plot(U[site, :], label=r"site_%d_iter_%d" % (site, i))
            plt.xlabel("Time")
            plt.ylabel(r"$U$")
            plt.title(r"Trajectory of U")
            legend = plt.legend(
                bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0
            )
            fig.tight_layout()
            plt.subplots_adjust(right=0.7)
            fig.savefig(
                path / f"U_site_{site}_iter_{i}.png",
                bbox_extra_artists=(legend,),
                bbox_inches="tight",
                dpi=100,
            )
            plt.close()


def coef_prior_density(coef_samples, plots_dir, K=8):
    theta_coef_samples = coef_samples[:, :K]
    gamma_coef_samples = coef_samples[:, K:]

    for i in range(theta_coef_samples.shape[1]):
        # Make paths
        path = plots_dir / "theta_coef_prior"
        if not os.path.exists(path):
            os.makedirs(path)

        fig, axes = plt.subplots(1, 1, figsize=(8, 6))
        plt.hist(theta_coef_samples[:, i], bins=30, alpha=0.7)
        plt.ylabel(r"$Frequency$")
        plt.title(r"Prior density of $\beta^\theta_%d$" % i)
        fig.tight_layout()
        plt.subplots_adjust(right=0.7)
        fig.savefig(
            path / f"theta_coef_{i}.png",
            bbox_inches="tight",
            dpi=100,
        )
        plt.close()

    for i in range(gamma_coef_samples.shape[1]):
        # Make paths
        path = plots_dir / "gamma_coef_prior"
        if not os.path.exists(path):
            os.makedirs(path)

        fig, axes = plt.subplots(1, 1, figsize=(8, 6))
        plt.hist(gamma_coef_samples[:, i], bins=30, alpha=0.7)
        plt.ylabel(r"$Frequency$")
        plt.title(r"Prior density of $\beta^\gamma_%d$" % i)
        fig.tight_layout()
        plt.subplots_adjust(right=0.7)
        fig.savefig(
            path / f"gamma_coef_{i}.png",
            bbox_inches="tight",
            dpi=100,
        )
        plt.close()


def coef_posterior_density(coef_samples, plots_dir, K=8):
    theta_coef_samples = coef_samples[:, :K]
    gamma_coef_samples = coef_samples[:, K:]

    for i in range(theta_coef_samples.shape[1]):
        # Make paths
        path = plots_dir / "theta_coef_posterior"
        if not os.path.exists(path):
            os.makedirs(path)

        fig, axes = plt.subplots(1, 1, figsize=(8, 6))
        plt.hist(theta_coef_samples[:, i], bins=30, alpha=0.7)
        plt.ylabel(r"$Frequency$")
        plt.title(r"Posterior density of $\beta^\theta_%d$" % i)
        fig.tight_layout()
        plt.subplots_adjust(right=0.7)
        fig.savefig(
            path / f"theta_coef_{i}.png",
            bbox_inches="tight",
            dpi=100,
        )
        plt.close()

    for i in range(gamma_coef_samples.shape[1]):
        # Make paths
        path = plots_dir / "gamma_coef_posterior"
        if not os.path.exists(path):
            os.makedirs(path)

        fig, axes = plt.subplots(1, 1, figsize=(8, 6))
        plt.hist(gamma_coef_samples[:, i], bins=30, alpha=0.7)
        plt.ylabel(r"$Frequency$")
        plt.title(r"Posterior density of $\beta^\gamma_%d$" % i)
        fig.tight_layout()
        plt.subplots_adjust(right=0.7)
        fig.savefig(
            path / f"gamma_coef_{i}.png",
            bbox_inches="tight",
            dpi=100,
        )
        plt.close()
