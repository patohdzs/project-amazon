import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.integrate import quad
from scipy.stats import truncnorm


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
            plt.xlabel("Iteration")
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


def Ua_trajectory(results: dict, plots_dir: Path) -> None:
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


def Um_trajectory(results: dict, plots_dir: Path) -> None:
    for site in range(results["size"]):
        path = plots_dir / "Um_trajectory" / f"site_{site}"
        if not os.path.exists(path):
            os.makedirs(path)

        for i, Um in enumerate(results["sol_val_Um_tracker"]):
            fig, axes = plt.subplots(1, 1, figsize=(8, 6))
            plt.plot(Um[site, :], label=r"site_%d_iter_%d" % (site, i))
            plt.xlabel("Iteration")
            plt.ylabel(r"$Um$")
            plt.title(r"Trajectory of Um")
            legend = plt.legend(
                bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0
            )
            fig.tight_layout()
            plt.subplots_adjust(right=0.7)
            fig.savefig(
                path / f"Um_site_{site}_iter_{i}.png",
                bbox_extra_artists=(legend,),
                bbox_inches="tight",
                dpi=100,
            )
            plt.close()


def Up_trajectory(results: dict, plots_dir: Path) -> None:
    for site in range(results["size"]):
        path = plots_dir / "Up_trajectory" / f"site_{site}"
        if not os.path.exists(path):
            os.makedirs(path)

        for i, Up in enumerate(results["sol_val_Up_tracker"]):
            fig, axes = plt.subplots(1, 1, figsize=(8, 6))
            plt.plot(Up[site, :], label=r"site_%d_iter_%d" % (site, i))
            plt.xlabel("Iteration")
            plt.ylabel(r"$Up$")
            plt.title(r"Trajectory of Up")
            legend = plt.legend(
                bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0
            )
            fig.tight_layout()
            plt.subplots_adjust(right=0.7)
            fig.savefig(
                path / f"Up_site_{site}_iter_{i}.png",
                bbox_extra_artists=(legend,),
                bbox_inches="tight",
                dpi=100,
            )
            plt.close()


def traceplot_Up(results, plotdir):
    for j in range(results["size"]):
        for i in range(len(results["sol_val_Up_tracker"])):
            i = len(results["sol_val_Up_tracker"]) - 1
            fig, axes = plt.subplots(1, 1, figsize=(8, 6))
            plt.plot(
                results["sol_val_Up_tracker"][i][j, :],
                label=r"$Up_{site_%d, iter_%d}$" % (j + 1, i),
            )
            plt.xlabel("Iteration")
            plt.ylabel(r"$Up$")
            plt.title(r"Trace Plot of Up")
            legend = plt.legend(
                bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0
            )
            fig.tight_layout()
            plt.subplots_adjust(right=0.7)
            fig.savefig(
                plotdir / "Up_site_%d_iter_%d.png" % (j + 1, i),
                bbox_extra_artists=(legend,),
                bbox_inches="tight",
                dpi=100,
            )
    plt.close()


def dist_theta(results, plotdir):
    size = results["size"]
    for j in range(size):
        len(
            results["collected_ensembles"][len(results["collected_ensembles"]) - 1][
                :, j
            ]
        )
        mu = theta[j]  # Assuming theta is your mean array
        sigma = thetaSD[j]  # Assuming thetaSD is your standard deviation array
        theta_samples = np.random.normal(mu, sigma, 50000)

        mu = gamma[j]  # Assuming gamma is your mean array
        sigma = gammaSD[j]  # Assuming gammaSD is your standard deviation array
        gamma_samples = np.random.normal(mu, sigma, 50000)

        # Find min and max across all iterations for the same site
        min_theta = np.min(
            [
                np.min(results["collected_ensembles"][i][:, j])
                for i in range(len(results["collected_ensembles"]))
            ]
            + [np.min(theta_samples)]
        )
        max_theta = np.max(
            [
                np.max(results["collected_ensembles"][i][:, j])
                for i in range(len(results["collected_ensembles"]))
            ]
            + [np.max(theta_samples)]
        )

        min_gamma = np.min(
            [
                np.min(results["collected_ensembles"][i][:, j + size])
                for i in range(len(results["collected_ensembles"]))
            ]
            + [np.min(gamma_samples)]
        )
        max_gamma = np.max(
            [
                np.max(results["collected_ensembles"][i][:, j + size])
                for i in range(len(results["collected_ensembles"]))
            ]
            + [np.max(gamma_samples)]
        )

        # for i in range(len(results['collected_ensembles'])):
        i = len(results["sol_val_Z_tracker"]) - 1

        # i = 0

        # For theta parameters
        fig, ax1 = plt.subplots(1, 1, figsize=(8, 6))
        sns.histplot(
            theta_samples[theta_samples > 0],
            bins=100,
            label="Unadjusted",
            kde=False,
            color="blue",
            ax=ax1,
        )
        ax1.set_ylabel("Count (Unadjusted)")
        ax2 = ax1.twinx()
        sns.histplot(
            results["collected_ensembles"][i][:, j][
                results["collected_ensembles"][i][:, j] > 0
            ],
            bins=100,
            label="Adjusted",
            kde=False,
            color="red",
            ax=ax2,
        )
        ax2.set_ylabel("Count (Adjusted)")
        plt.xlim([min_theta, max_theta])  # Set x-axis limits
        plt.xlabel(r"$\theta_{site_%d}$" % (j + 1))
        # plt.title(r"Distribution of $\theta_{site_%d}$ for iteration %d"%(j+1, i))
        plt.title(
            r"Distribution of $\theta_{site_%d}$ reflection boundary methods" % (j + 1)
        )
        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, loc=0)
        ax1.grid(False)
        ax2.grid(False)
        fig.savefig(plotdir + "theta_site_%d_iter_%d_hist.png" % (j + 1, i), dpi=100)
        plt.close()

        fig, ax1 = plt.subplots(1, 1, figsize=(8, 6))
        sns.histplot(
            gamma_samples[gamma_samples > 0],
            bins=100,
            label="Unadjusted",
            kde=False,
            color="blue",
            ax=ax1,
        )
        ax1.set_ylabel("Count (Unadjusted)")
        ax2 = ax1.twinx()
        sns.histplot(
            results["collected_ensembles"][i][:, j + size][
                results["collected_ensembles"][i][:, j + size] > 0
            ],
            bins=100,
            label="Adjusted",
            kde=False,
            color="red",
            ax=ax2,
        )
        ax2.set_ylabel("Count (Adjusted)")
        plt.xlim([min_gamma, max_gamma])  # Set x-axis limits
        plt.xlabel(r"$\gamma_{site_%d}$" % (j + 1))
        plt.title(
            r"Distribution of $\gamma_{site_%d}$ for reflection boundary methods"
            % (j + 1)
        )
        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, loc=0)
        ax1.grid(False)
        ax2.grid(False)
        fig.savefig(plotdir + "gamma_site_%d_iter_%d_hist.png" % (j + 1, i), dpi=100)
        plt.close()

        fig, axes = plt.subplots(1, 1, figsize=(8, 6))
        # sns.kdeplot(
        #     theta_samples, label="Unadjusted", shade=True, clip=(0, None),
        #     color="blue"
        # )
        a, b = (0 - np.mean(theta_samples)) / np.std(theta_samples), np.inf
        unadjusted_dist = truncnorm(
            a, b, loc=np.mean(theta_samples), scale=np.std(theta_samples)
        )
        x = np.linspace(min_theta, max_theta, 1000)
        axes.plot(
            x, unadjusted_dist.pdf(x), label="Unadjusted (truncnorm)", color="blue"
        )  # plot scipy result as dashed
        axes.fill_between(
            x, unadjusted_dist.pdf(x), color="blue", alpha=0.3
        )  # Shaded area for unadjusted
        sns.kdeplot(
            results["collected_ensembles"][i][:, j],
            label="Adjusted",
            shade=True,
            clip=(0, None),
            color="red",
        )
        plt.xlim([min_theta, max_theta])  # Set x-axis limits
        plt.xlabel(r"$\theta_{site_%d}$" % (j + 1))
        plt.ylabel("Density")
        # plt.title(r"Density of $\theta_{site_%d}$ for iteration %d"%(j+1, i))
        plt.title(
            r"Density of $\theta_{site_%d}$ for reflection boundary methods" % (j + 1)
        )
        if j in [0, 1]:
            plt.ylim([0, 120])
        plt.legend()
        fig.savefig(plotdir + "theta_site_%d_iter_%d.png" % (j + 1, i), dpi=100)
        plt.close()

        fig, axes = plt.subplots(1, 1, figsize=(8, 6))
        # sns.kdeplot(
        #     gamma_samples, label="Unadjusted", shade=True, clip=(0, None),
        #     color="blue"
        # )
        a, b = (0 - np.mean(gamma_samples)) / np.std(gamma_samples), np.inf
        truncnorm(a, b, loc=np.mean(gamma_samples), scale=np.std(gamma_samples))
        x = np.linspace(min_theta, max_theta, 1000)
        axes.plot(
            x, unadjusted_dist.pdf(x), label="Unadjusted (truncnorm)", color="blue"
        )  # plot scipy result as dashed
        axes.fill_between(
            x, unadjusted_dist.pdf(x), color="blue", alpha=0.3
        )  # Shaded area for unadjusted
        sns.kdeplot(
            results["collected_ensembles"][i][:, j + size],
            label="Adjusted",
            shade=True,
            clip=(0, None),
            color="red",
        )
        plt.xlim([min_gamma, max_gamma])  # Set x-axis limits
        plt.xlabel(r"$\gamma_{site_%d}$" % (j + 1))
        plt.ylabel("Density")
        # plt.title(r"Density of $\gamma_{site_%d}$ for iteration %d"%(j+1, i))
        plt.title(
            r"Density of $\gamma_{site_%d}$ for reflection boundary methods" % (j + 1)
        )
        plt.legend()
        fig.savefig(plotdir + "gamma_site_%d_iter_%d.png" % (j + 1, i), dpi=100)
        plt.close()

        fig, axes = plt.subplots(1, 1, figsize=(8, 6))
        # sns.kdeplot(
        #     theta_samples, label="Unadjusted", shade=True, clip=(0, None),
        #    color="blue"
        # )
        # sns.kdeplot(
        #     results["collected_ensembles"][i][:, j],
        #     label="Adjusted",
        #     shade=True,
        #     clip=(0, None),
        #     color="red",
        # )
        a, b = (0 - np.mean(theta_samples)) / np.std(theta_samples), np.inf
        unadjusted_dist = truncnorm(
            a, b, loc=np.mean(theta_samples), scale=np.std(theta_samples)
        )
        a, b = (0 - np.mean(results["collected_ensembles"][i][:, j])) / np.std(
            results["collected_ensembles"][i][:, j]
        ), np.inf
        adjusted_dist = truncnorm(
            a,
            b,
            loc=np.mean(results["collected_ensembles"][i][:, j]),
            scale=np.std(results["collected_ensembles"][i][:, j]),
        )
        x = np.linspace(min_theta, max_theta, 1000)
        axes.plot(
            x, unadjusted_dist.pdf(x), label="Unadjusted (truncnorm)", color="blue"
        )  # plot scipy result as dashed
        axes.fill_between(
            x, unadjusted_dist.pdf(x), color="blue", alpha=0.3
        )  # Shaded area for unadjusted
        axes.plot(
            x, adjusted_dist.pdf(x), label="Adjusted (truncnorm)", color="red"
        )  # plot scipy result as dashed
        axes.fill_between(x, adjusted_dist.pdf(x), color="red", alpha=0.3)
        unadjusted_integral, _ = quad(unadjusted_dist.pdf, min_theta, max_theta)
        adjusted_integral, _ = quad(adjusted_dist.pdf, min_theta, max_theta)
        axes.text(
            0.6,
            0.6,
            f"""
            Unadjusted integral: {unadjusted_integral:.2f}\n
            Adjusted integral: {adjusted_integral:.2f}
            """,
            transform=axes.transAxes,
        )

        plt.xlim([min_theta, max_theta])  # Set x-axis limits
        plt.xlabel(r"$\theta_{site_%d}$" % (j + 1))
        plt.ylabel("Density")
        # plt.title(r"Density of $\theta_{site_%d}$ for iteration %d"%(j+1, i))
        plt.title(
            r"Density of $\theta_{site_%d}$ for reflection boundary methods" % (j + 1)
        )
        if j in [0, 1]:
            plt.ylim([0, 120])
        plt.legend()
        fig.savefig(
            plotdir + "theta_site_%d_iter_%d_truncnormal.png" % (j + 1, i), dpi=100
        )
        plt.close()

        fig, axes = plt.subplots(1, 1, figsize=(8, 6))
        # sns.kdeplot(
        #     gamma_samples, label="Unadjusted", shade=True, clip=(0, None),
        #     color="blue"
        # )
        # sns.kdeplot(
        #     results["collected_ensembles"][i][:, j + size],
        #     label="Adjusted",
        #     shade=True,
        #     clip=(0, None),
        #     color="red",
        # )
        a, b = (0 - np.mean(gamma_samples)) / np.std(gamma_samples), np.inf
        truncnorm(a, b, loc=np.mean(gamma_samples), scale=np.std(gamma_samples))
        a, b = (0 - np.mean(results["collected_ensembles"][i][:, j + size])) / np.std(
            results["collected_ensembles"][i][:, j + size]
        ), np.inf
        truncnorm(
            a,
            b,
            loc=np.mean(results["collected_ensembles"][i][:, j + size]),
            scale=np.std(results["collected_ensembles"][i][:, j + size]),
        )
        x = np.linspace(min_gamma, max_gamma, 1000)
        axes.fill_between(
            x, unadjusted_dist.pdf(x), alpha=0.5, color="blue"
        )  # plot scipy result as filled area
        axes.plot(
            x, unadjusted_dist.pdf(x), label="Unadjusted (truncnorm)", color="blue"
        )  # plot scipy result as solid line
        axes.fill_between(
            x, adjusted_dist.pdf(x), alpha=0.5, color="red"
        )  # plot scipy result as filled area
        axes.plot(
            x, adjusted_dist.pdf(x), label="Adjusted (truncnorm)", color="red"
        )  # plot scipy result as solid line
        unadjusted_integral, _ = quad(unadjusted_dist.pdf, min_gamma, max_gamma)
        adjusted_integral, _ = quad(adjusted_dist.pdf, min_gamma, max_gamma)
        axes.text(
            0.6,
            0.6,
            f"""
            Unadjusted integral: {unadjusted_integral:.2f}\n
            Adjusted integral: {adjusted_integral:.2f}
            """,
            transform=axes.transAxes,
        )

        plt.xlim([min_gamma, max_gamma])  # Set x-axis limits
        plt.xlabel(r"$\gamma_{site_%d}$" % (j + 1))
        plt.ylabel("Density")
        # plt.title(r"Density of $\gamma_{site_%d}$ for iteration %d"%(j+1, i))
        plt.title(
            r"Density of $\gamma_{site_%d}$ for reflection boundary methods" % (j + 1)
        )
        plt.legend()
        fig.savefig(
            plotdir + "gamma_site_%d_iter_%d_truncnormal.png" % (j + 1, i), dpi=100
        )
        plt.close()
