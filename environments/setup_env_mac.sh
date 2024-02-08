#!/usr/bin/env bash
set -e
conda create --yes --name ahtvs psycopg2 matplotlib pandas jupyter python=3.6 #issues with newer versions of rdkit and python currently on macs
conda install conda-build
conda activate ahtvs
conda install --yes -c rdkit rdkit=2019 #issues with newer versions of rdkit currently on macs
conda install --yes -c conda-forge --no-channel-priority django=2.0 django-debug-toolbar django-extensions django-filter django-guardian djangorestframework=3.9.3 # newer versions of django causing issues
conda install -c openbabel openbabel
pip install --upgrade pip
pip install -r ahtvs-docker-dev/requirements.txt
pip install -e ../djangochem
conda develop ../
#conda install -c anaconda conda-build=2.1.10 #Issues with newer versions of conda-build looking for the wrong path, no longer needed Jan 7 2020
#conda develop ../
#conda develop ../djangochem
#pip install -e ../
#pip install -e ../djangochem
