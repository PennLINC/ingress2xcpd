# Ingress2QSIRecon
A tool for ingressing outputs from other processing pipelines (e.g., HCP, UKB) for use in QSIRecon.

## Overview
This tool can be used to import data from bespoke and widely-used DWI preprocessing pipelines (e.g., Human Connectome Project and UK Biobank) into QSIRecon for post-processing.

## Installation
This project can be installed via PyPi. In your virtual environment run
```
pip install Ingress2QSIRecon
```
For the master/development branch, you can clone the github repo, and run from within it:
```
pip install -e .
```

## Usage
Assuming you have the data to ingress locally availble, you can run the following:
```
ingress2qsirecon \
    /PATH/TO/INPUT/DATA \ # E.g., path to HCP1200 directory
    /PATH/TO/OUTPUT \ # BIDS directory with ingressed data will be made here
    PIPELINE_NAME \ # Currently support "hcpya" and "ukb" \
    -w /PATH/TO/WORKING/DIR \ # Where intermediate files and workflow information will be stored
    --participant-label 100307 # Name(s) of folder within /PATH/TO/INPUT/DATA. If not specified, will process everyone
```