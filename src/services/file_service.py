import os
from pathlib import Path


def output_dir_path(**kwargs):
    path = get_path("output")
    return _nested_dir_path(path, **kwargs)


def logs_dir_path(**kwargs):
    path = get_path("logs")
    return _nested_dir_path(path, **kwargs)


def plots_dir_path(**kwargs):
    path = get_path("plots")
    return _nested_dir_path(path, **kwargs)


def get_path(*args: str) -> Path:
    return _project_root().joinpath(*args)


def _nested_dir_path(path, **kwargs):
    for key, value in kwargs.items():
        # Make subdir name
        subdir_name = f"{key}_{value}"

        # Join subdir to rest of the path
        path = path.joinpath(path, subdir_name)

        # Create subdir if non-existent
        if not os.path.exists(path):
            os.makedirs(path)

    return path


def _project_root() -> Path:
    return Path(__file__).parents[2]
