#!/usr/bin/env bash
set -e
conda create --yes --name ahtvs psycopg2 matplotlib pandas jupyter python=3
conda activate ahtvs
conda install --yes -c rdkit rdkit #issues with newer versions of rdkit currently on macs
conda install --yes -c conda-forge --no-chan-pri django django-debug-toolbar django-extensions django-filter django-guardian djangorestframework==3.9.4 # problems with six djangorestframework past this one
pip install --upgrade pip
pip install -r ahtvs-docker-dev/requirements.txt
conda install conda-build=2.1.10 #Issues with newer versions of conda-build looking for the wrong path
conda develop ../djangochem
conda develop ../
#pip install -e ../djangochem
#pip install -e ../
