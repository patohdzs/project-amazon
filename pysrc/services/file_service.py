import os
from pathlib import Path


def output_dir_path(**kwargs) -> Path:
    """Get path to the output directory.

    Returns:
        Path: Path to the output directory.
    """
    path = get_path("output")
    return _nested_dir_path(path, **kwargs)


def logs_dir_path(**kwargs) -> Path:
    """Get path to the logs directory.

    Returns:
        Path: Path to the logs directory.
    """
    path = get_path("logs")
    return _nested_dir_path(path, **kwargs)


def plots_dir_path(**kwargs) -> Path:
    """Get path to the plots directory.

    Returns:
        Path: Path to the plots directory.
    """
    path = get_path("plots")
    return _nested_dir_path(path, **kwargs)


def stan_model_path(model_name: str) -> Path:
    return get_path("stan_models", model_name)


def get_path(*args: str) -> Path:
    return _project_root().joinpath(*args)


def _nested_dir_path(path, **kwargs) -> Path:
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
