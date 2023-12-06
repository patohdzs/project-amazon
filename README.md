# Project Amazon

## Requirements
- Python >= 3.9
- Project data with the following directory structure:

```
.
└── data
    ├── calibration
    │   ├── farm_gate_price.xlsx
    │   ├── ipeadata[21-08-2023-01-28].xls
    │   └── prepData
    │       ├── muniTheta_prepData.Rdata
    │       ├── muniTheta_prepData_gamma.Rdata
    │       └── seriesPriceCattle_prepData.Rdata
    └── hmc
        ├── data_gamma.geojson
        ├── data_theta.geojson
        ├── hmc_10SitesModel.csv
        ├── hmc_24SitesModel.csv
        ├── hmc_40SitesModel.csv
        ├── id_10.geojson
        ├── id_24.geojson
        └── id_40.geojson
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

3. Install Stan
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
