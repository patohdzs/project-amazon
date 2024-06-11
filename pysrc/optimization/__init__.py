import pandas as pd

from pysrc.services.file_service import get_path


def save_trajectories(trajectories):
    for var, traj in trajectories.items():
        path = get_path("data", "trajectories") / f"amazon_data_{var.lower()}.dat"
        pd.DataFrame(traj).to_csv(path, sep="\t", index=False)
