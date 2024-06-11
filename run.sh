#!/bin/bash

# Suggested usage function
usage() {
    echo "Usage: $0 [-spcdra]" 1>&2
    exit 1
}

setup_flag='false'
processing_flag='false'  # data/raw -> data/clean
calibration_flag='false' # data/clean -> data/calibration
det_model_flag='false'   # data/calibration -> data/trajectories/det
mpc_model_flag='false'   # data/calibration -> data/trajectories/mpc
hmc_model_flag='false'   # data/calibration -> data/trajectories/hmc
analysis_flag='false'    # data/trajectories -> output/

# Check which tasks to run
while getopts 'spcdmha' flag; do
    case "${flag}" in
    s) setup_flag='true' ;;
    p) processing_flag='true' ;;
    c) calibration_flag='true' ;;
    d) det_model_flag='true' ;;
    m) mpc_model_flag='true' ;;
    h) hmc_model_flag='true' ;;
    a) analysis_flag='true' ;;
    *) usage ;;
    esac
done

if [ "$setup_flag" = "true" ]; then
    # Check if Python is installed
    if ! command -v python &>/dev/null; then
        echo "Python is not installed. Please install Python first."
        exit 1
    fi

    # Check if Python is installed
    if ! command -v Rscript &>/dev/null; then
        echo "R is not installed. Please install R first."
        exit 1
    fi

    # Check if one of Gurobi or GAMs is installed

    # Check if venv is installed
    if ! command -v python3 -m venv &>/dev/null; then
        echo "venv is not installed. Install by running 'python3 -m pip install --user virtualenv'."
        exit 1
    fi

    # Create virtual environment
    echo Creating python virtual environment...
    python3 -m venv .venv

    # Activate virtual environment
    source .venv/bin/activate
    echo Done!

    # Install python dependencies
    echo Installing python dependencies...
    python -m pip install -e '.[all]'
    echo Done!

fi

if [ "$processing_flag" = "true" ]; then
    echo "Cleaning and processing data..."
    Rscript rsrc/raw2clean/_masterfile_raw2clean.R
    echo "Done!"
fi

if [ "$calibration_flag" = "true" ]; then
    echo "Calibrating model parameters..."
    Rscript rsrc/calibration/_masterfile_prep.R
    Rscript rsrc/calibration/_masterfile_calibration.R
    echo "Done!"
fi

if [ "$det_model_flag" = "true" ]; then
    # Running deterministic model script
    echo Running deterministic model...
    python3 scripts/run_det_model.py
    echo Done!
fi

if [ "$mpc_model_flag" = "true" ]; then
    # Running price risk model script
    echo Running deterministic model...
    python3 scripts/run_mpc_model.py
    echo Done!
fi

if [ "$hmc_model_flag" = "true" ]; then
    # Running ambiguity model script
    echo Running deterministic model...
    python3 scripts/run_hmc_model.py
    echo Done!
fi
