# Project Amazon

## Requirements
- Python >= 3.9, <3.12
- Gurobi >= 10.0.3
## Data Requirements
To replicate, make sure to download the raw data into the directory structure below:
```
.
└── data
    └── raw
        ├── esa
        ├── fgv
        ├── global_forest_watch
        ├── ibge
        ├── ipea
        ├── mapbiomas
        ├── seabpr
        ├── seeg
        ├── worldbank
        └── worldclim
```

## Installation

0. Clone git repository and move into `project-amazon/`
1. Create and activate a new virtual environment
```
python -m venv .venv
source .venv/bin/activate
```
2. Install python dependencies
```
python -m pip install -e '.[all]'
```

3. Install CmdStan
```
install_cmdstan --overwrite
```

4. Install pre-commit hooks (required for contributors)
```
pre-commit install
```


## Contributing
0. Open a new git branch
```
git checkout -b <new branch name>
```
1. Create new code changes

2. Stage changed files
```
git add <names of changed files>
```

3. Commit changed files
```
git commit
```

4. After several commits, push commits to remote (if it is the first time pushing this branch use the `--set-upstream` flag)
```
git push
```

5. Submit a pull request on GitHub
