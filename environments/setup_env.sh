#!/usr/bin/env bash
set -e
conda create --yes --name ahtvs python=3
source activate ahtvs
conda install --yes psycopg2 matplotlib pandas jupyter
conda install --yes -c rdkit rdkit
conda install --yes -c conda-forge --no-chan-pri django django-debug-toolbar django-extensions django-filter django-guardian djangorestframework
pip install --upgrade pip
pip install -r ahtvs-docker-dev/requirements.txt
pip install -e ../../djangochem


