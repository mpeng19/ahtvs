#!/bin/bash

# custom location to install the a2g2 environments
# full paths are recommended
# if you use '~', don't put quotes around it
ENV_FOLDER=

# create environment folder
if [ -z "$ENV_FOLDER" ]; then
    ENV_FOLDER=~/a2g2-singularity
fi
mkdir -p $ENV_FOLDER

# obtaining environment variables
SRCDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DCHEM_PATH=../../djangochem
if [ -d "$DCHEM_PATH" ]; then
    DCHEM_PATH="$(cd $DCHEM_PATH; pwd -P)"
else
    DCHEM_PATH="/your/path/to/djangochem"
fi
cd $ENV_FOLDER

# pull postgres database
if [ -f postgres.img ]; then
    echo 'Found postgres.img'
else
    echo 'Download postgres.img'
    singularity pull --name postgres.img shub://verysure/postgres-alpine
fi

# pull rdkit
if [ -f rdkit-django.img ]; then
    echo 'Found rdkit-django.img'
else
    echo 'Download rdkit-django.img'
    singularity pull --name rdkit-django.img shub://verysure/rdkit-django
fi

# create settings file
if [ ! -f settings ]; then
    cp $SRCDIR/settings ./
    echo "export PYTHONPATH=$DCHEM_PATH" >> settings
    echo "Created settings:\n\tRemember to rename database, password and a2g2 djangochem path (PYTHONPATH)"
fi

# create local host pgdata folder
if [ ! -d ./pgdata ]; then
    mkdir -p ./pgdata
fi

# copy run script
if [ ! -f a2g2 ]; then
    cp $SRCDIR/a2g2 ./
    chmod +x ./a2g2
    echo "Copied executable a2g2, use $ENV_FOLDER/a2g2 for starting database, open shell and more."
    echo "Include alias in bash for easier usage: alias a2g2=$ENV_FOLDER/a2g2"
fi

