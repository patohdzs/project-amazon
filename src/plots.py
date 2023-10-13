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


def traceplot_sampling_time(results: dict, plots_dir: Path) -> None:
    # Make plot fig and axis
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    plt.plot(results["sampling_time_tracker"], label=r"Sampling Time")
    plt.xlabel("Iteration")
    plt.ylabel("Sampling Time (s)")
    plt.title(r"Trace Plot of Sampling Time")
    legend = plt.legend(bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0)
    fig.tight_layout()
    plt.subplots_adjust(right=0.7)
    fig.savefig(
        plots_dir / "sampling_time.png",
        bbox_extra_artists=(legend,),
        bbox_inches="tight",
        dpi=100,
    )
    plt.close()


def traceplot_params_pct_error(results: dict, plots_dir: Path) -> None:
    num_sites = results["num_sites"]
    # Find pct change
    pct_change = (
        abs(np.diff(results["uncertain_vals_tracker"], axis=0))
        / results["uncertain_vals_tracker"][:-1, :]
    )

    # Create mpl base
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    for i in range(num_sites):
        axes[0].plot(pct_change[:, i], label=f"Theta {i+1}")

    axes[0].set_title("Theta Proportional Error")
    for i in range(num_sites, pct_change.shape[1]):
        axes[1].plot(pct_change[:, i], label=f"Gamma {i+1-num_sites}")

    axes[1].set_title("Gamma Proportional Error")
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


def traceplot_gammas(results: dict, plots_dir: Path) -> None:
    num_sites = results["num_sites"]
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    for j in range(num_sites):
        plt.plot(
            results["uncertain_vals_tracker"][:, j + num_sites],
            label=r"$\gamma_{%d}$" % (j + 1),
        )
    plt.xlabel("Iteration")
    plt.ylabel(r"$\gamma$")
    plt.title(r"Trace Plot of $\gamma$")
    legend = plt.legend(bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0)
    fig.tight_layout()
    plt.subplots_adjust(right=0.7)
    fig.savefig(
        plots_dir / "gamma.png",
        bbox_extra_artists=(legend,),
        bbox_inches="tight",
        dpi=100,
    )
    plt.close()


def traceplot_thetas(results: dict, plots_dir: Path) -> None:
    num_sites = results["num_sites"]
    fig, axes = plt.subplots(1, 1, figsize=(8, 6))
    for j in range(num_sites):
        plt.plot(
            results["uncertain_vals_tracker"][:, j],
            label=r"$\theta_{%d}$" % (j + 1),
        )
    plt.xlabel("Iteration")
    plt.ylabel(r"$\theta$")
    plt.title(r"Trace Plot of $\theta$")
    legend = plt.legend(bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0)
    fig.tight_layout()
    plt.subplots_adjust(right=0.7)
    fig.savefig(
        plots_dir / "theta.png",
        bbox_extra_artists=(legend,),
        bbox_inches="tight",
        dpi=100,
    )
    plt.close()


def Z_trajectory(results: dict, plots_dir: Path) -> None:
    for site in range(results["num_sites"]):
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
    for site in range(results["num_sites"]):
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
    for site in range(results["num_sites"]):
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
    for site in range(results["num_sites"]):
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


def prior_density(theta_samples, gamma_samples, plots_dir, num_sites=10):
    for i in range(theta_samples.shape[1]):
        # Make paths
        path = plots_dir / "theta_prior"
        if not os.path.exists(path):
            os.makedirs(path)

        fig, axes = plt.subplots(1, 1, figsize=(8, 6))
        plt.hist(theta_samples[:, i], density=True, bins=60, alpha=0.7)
        plt.ylabel(r"$Frequency$")
        plt.title(r"Prior density of $\theta_%d$" % i)
        fig.tight_layout()
        plt.subplots_adjust(right=0.7)
        fig.savefig(
            path / f"theta_{i}.png",
            bbox_inches="tight",
            dpi=100,
        )
        plt.close()

    for i in range(gamma_samples.shape[1]):
        # Make paths
        path = plots_dir / "gamma_prior"
        if not os.path.exists(path):
            os.makedirs(path)

        fig, axes = plt.subplots(1, 1, figsize=(8, 6))
        plt.hist(gamma_samples[:, i], density=True, bins=60, alpha=0.7)
        plt.ylabel(r"$Frequency$")
        plt.title(r"Prior density of $\gamma_%d$" % i)
        fig.tight_layout()
        plt.subplots_adjust(right=0.7)
        fig.savefig(
            path / f"gamma_{i}.png",
            bbox_inches="tight",
            dpi=100,
        )
        plt.close()


def posterior_density(samples, plots_dir, num_sites=10):
    theta_samples = samples[:, :num_sites]
    gamma_samples = samples[:, num_sites:]

    for i in range(theta_samples.shape[1]):
        # Make paths
        path = plots_dir / "theta_posterior"
        if not os.path.exists(path):
            os.makedirs(path)

        fig, axes = plt.subplots(1, 1, figsize=(8, 6))
        plt.hist(theta_samples[:, i], density=True, bins=30, alpha=0.7)
        plt.ylabel(r"$Frequency$")
        plt.title(r"Posterior density of $\theta_%d$" % i)
        fig.tight_layout()
        plt.subplots_adjust(right=0.7)
        fig.savefig(
            path / f"theta_{i}.png",
            bbox_inches="tight",
            dpi=100,
        )
        plt.close()

    for i in range(gamma_samples.shape[1]):
        # Make paths
        path = plots_dir / "gamma_posterior"
        if not os.path.exists(path):
            os.makedirs(path)

        fig, axes = plt.subplots(1, 1, figsize=(8, 6))
        plt.hist(gamma_samples[:, i], density=True, bins=30, alpha=0.7)
        plt.ylabel(r"$Frequency$")
        plt.title(r"Posterior density of $\gamma_%d$" % i)
        fig.tight_layout()
        plt.subplots_adjust(right=0.7)
        fig.savefig(
            path / f"gamma_{i}.png",
            bbox_inches="tight",
            dpi=100,
        )
        plt.close()


def overlap_prior_posterior(prior_samples, post_samples, plots_dir, num_sites=10):
    theta_prior_samples = prior_samples[:, :num_sites]
    gamma_prior_samples = prior_samples[:, num_sites:]

    theta_post_samples = post_samples[:, :num_sites]
    gamma_post_samples = post_samples[:, num_sites:]

    for i in range(theta_prior_samples.shape[1]):
        # Make paths
        path = plots_dir / "theta_density"
        if not os.path.exists(path):
            os.makedirs(path)

        prior = theta_prior_samples[:, i]
        post = theta_post_samples[:, i]
        bins = np.histogram(np.hstack((prior, post)), bins=60)[1]

        fig, axes = plt.subplots(1, 1, figsize=(8, 6))
        plt.hist(prior, density=True, bins=bins, alpha=0.7)
        plt.hist(post, density=True, bins=bins, alpha=0.7, color="red")
        plt.ylabel(r"$Frequency$")
        plt.title(r"Density of $\theta_%d$" % i)
        fig.tight_layout()
        plt.subplots_adjust(right=0.7)
        fig.savefig(
            path / f"theta_{i}.png",
            bbox_inches="tight",
            dpi=100,
        )
        plt.close()

    for i in range(gamma_prior_samples.shape[1]):
        # Make paths
        path = plots_dir / "gamma_density"
        if not os.path.exists(path):
            os.makedirs(path)

        prior = gamma_prior_samples[:, i]
        post = gamma_post_samples[:, i]
        bins = np.histogram(np.hstack((prior, post)), bins=60)[1]

        fig, axes = plt.subplots(1, 1, figsize=(8, 6))
        plt.hist(prior, density=True, bins=bins, alpha=0.7)
        plt.hist(post, density=True, bins=bins, alpha=0.7, color="red")
        plt.ylabel(r"$Frequency$")
        plt.title(r"Density of $\gamma_%d$" % i)
        fig.tight_layout()
        plt.subplots_adjust(right=0.7)
        fig.savefig(
            path / f"gamma_{i}.png",
            bbox_inches="tight",
            dpi=100,
        )
        plt.close()
