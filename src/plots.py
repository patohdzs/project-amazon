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


def traceplot_params_pct_error(results: dict, plots_dir: Path) -> None:
    size = results["size"]
    diff = (
        abs(np.diff(results["uncertain_vals_tracker"], axis=0))
        / results["uncertain_vals_tracker"][:-1, :]
    )
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    for i in range(size):
        axes[0].plot(diff[:, i], label=f"Theta {i+1}")
    axes[0].set_title("Theta Proportional Error")
    for i in range(size, diff.shape[1]):
        axes[1].plot(diff[:, i], label=f"Gamma {i+1-size}")
    axes[1].set_title("Gamma Proportional Error")
    y_min = np.min(diff)
    y_max = np.max(diff)
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


def traceplot_gamma(results: dict, plots_dir: Path) -> None:
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    for j in range(results["size"], results["size"] * 2):
        plt.plot(
            results["uncertain_vals_tracker"][:, j], label=r"$\gamma_{%d}$" % (j + 1)
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


def traceplot_theta(results: dict, plots_dir: Path) -> None:
    fig, axes = plt.subplots(1, 1, figsize=(8, 6))
    for j in range(results["size"]):
        plt.plot(
            results["uncertain_vals_tracker"][:, j], label=r"$\theta_{%d}$" % (j + 1)
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


def traceplot_state(results: dict, plots_dir: Path) -> None:
    fig, axes = plt.subplots(1, 1, figsize=(8, 6))
    for j in range(results["size"] + 2):
        plt.plot(results["sol_val_X_tracker"][:, j], label=r"$X_{%d}$" % (j + 1))
    plt.xlabel("Iteration")
    plt.ylabel(r"$X$")
    plt.title(r"Trace Plot of X")
    legend = plt.legend(bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0)
    fig.tight_layout()
    plt.subplots_adjust(right=0.7)
    fig.savefig(
        plots_dir / "X.png", bbox_extra_artists=(legend,), bbox_inches="tight", dpi=100
    )
    plt.close()


def traceplot_adjustments(results: dict, plots_dir: Path) -> None:
    fig, axes = plt.subplots(1, 1, figsize=(8, 6))
    for j in range(results["size"] + 2):
        plt.plot(results["sol_val_Ua_tracker"][:, j], label=r"$Ua_{%d}$" % (j + 1))
    plt.xlabel("Iteration")
    plt.ylabel(r"$Ua$")
    plt.title(r"Trace Plot of Ua")
    legend = plt.legend(bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0)
    fig.tight_layout()
    plt.subplots_adjust(right=0.7)
    fig.savefig(
        plots_dir / "Ua.png", bbox_extra_artists=(legend,), bbox_inches="tight", dpi=100
    )
    plt.close()


def traceplot_Um(results: dict, plots_dir: Path) -> None:
    size = results["size"]
    for j in range(size):
        for i in range(len(results["sol_val_Um_tracker"])):
            i = len(results["sol_val_Up_tracker"]) - 1
            fig, axes = plt.subplots(1, 1, figsize=(8, 6))
            plt.plot(
                results["sol_val_Um_tracker"][i][j, :],
                label=r"site_%d_iter_%d" % (j + 1, i + 1),
            )
            plt.xlabel("Iteration")
            plt.ylabel(r"$Um$")
            plt.title(r"Trace Plot of Um for site_%d_iter_%d" % (j + 1, i + 1))
            legend = plt.legend(
                bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0
            )
            fig.tight_layout()
            plt.subplots_adjust(right=0.7)
            fig.savefig(
                plots_dir / "Um_site_%d_iter_%d.png" % (j + 1, i + 1),
                bbox_extra_artists=(legend,),
                bbox_inches="tight",
                dpi=100,
            )
            plt.close()
