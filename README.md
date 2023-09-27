# Project Amazon

## Requirements
- Python >= 3.9
- Project data with the following directory structure:

```
.
└── data
    └── hmc
        ├── calibration_10SitesModel.csv
        ├── calibration_24SitesModel.csv
        ├── calibration_40SitesModel.csv
        ├── data_gamma.geojson
        ├── data_theta.geojson
        ├── gamma_coe_ori.csv
        ├── gamma_vcov.csv
        ├── id_10.geojson
        ├── id_24.geojson
        ├── id_40.geojson
        ├── theta_coe_ori.csv
        └── theta_vcov.csv
```


## Installation

0. Clone git repository and move into `project-amazon/`
1. Create and activate a new virtual environment
```
python -m venv venv
source venv/bin/activate
```
2. Install dependencies
```
python -m pip install -e '.[all]'
```

3. Install pre-commit hooks (required for contributors)
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
