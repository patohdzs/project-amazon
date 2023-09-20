import os
from pathlib import Path


def get_output_dir(**kwargs):
    path = get_path("output")

    for key, value in kwargs.items():
        subdir_name = f"{key}_{value}"

        # Construct the full path of the subdirectory
        path = path.joinpath(path, subdir_name)

        # Create the subdirectory
        if not os.path.exists(path):
            os.makedirs(path)

    return path


def get_path(*args: str) -> Path:
    return _project_root().joinpath(*args)


def _project_root() -> Path:
    return Path(__file__).parents[2]
