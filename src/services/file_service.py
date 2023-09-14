from pathlib import Path


def get_path(*args: str) -> Path:
    return _project_root().joinpath(*args)


def _project_root() -> Path:
    return Path(__file__).parents[2]
