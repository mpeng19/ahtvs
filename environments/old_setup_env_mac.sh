#!/usr/bin/env bash
set -e
conda create --yes --name ahtvs psycopg2 matplotlib pandas jupyter python=3.5.3 #issues with newer versions of rdkit and python currently on macs
source activate ahtvs
conda install --yes -c rdkit rdkit=2017.03.03 #issues with newer versions of rdkit currently on macs
conda install --yes -c conda-forge --no-chan-pri django django-debug-toolbar django-extensions django-filter django-guardian djangorestframework
pip install --upgrade pip
pip install -r ahtvs-docker-dev/requirements.txt
pip install -e ../../djangochem
