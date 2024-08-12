# Project Amazon

## Requirements
- Python >= 3.9
- Gurobi >= 10.0.3 or GAMs >= (add version)
(Currently only works with Gurobi solver, please make sure install Gurobi before using the quick_script)

## Data directory structure
```
.
└── data
    └── hmc
        ├── muni_data_gamma.geojson
        ├── muni_data_theta.geojson
        ├── hmc_78SitesModel.csv
        ├── hmc_1043SitesModel.csv
        └── id_78.geojson
        └── id_1043.geojson
        └── site_78_data_gamma.geojson
        └── site_78_data_theta.geojson
        └── site_1043_data_gamma.geojson
        └── site_1043_data_theta.geojson
```

## Installation

0. Clone git repository and move into `project-amazon/`
1. Create and activate a new virtual environment
```
python -m venv venv
source venv/bin/activate
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
