BootStrap: docker
From: ubuntu:latest

%labels
MAINTAINER verysure

%files
files/initrc /
files/shell.sh /bin/

%post
apt-get -qq update --fix-missing 
apt-get install -yq wget bzip2 libxrender-dev libxext-dev libsm-dev
chmod a+x /bin/shell.sh

# install conda
wget --quiet https://repo.continuum.io/miniconda/Miniconda3-4.5.1-Linux-x86_64.sh -O /miniconda.sh
bash /miniconda.sh -b -p /opt/conda
rm /miniconda.sh

# install rdkit and django dependencies
. /opt/conda/etc/profile.d/conda.sh
conda activate base
conda install --yes conda
conda install --yes -c rdkit rdkit
conda install --yes psycopg2 pandas
conda install --yes -c conda-forge --no-chan-pri django django-debug-toolbar django-extensions django-filter django-guardian djangorestframework
pip install djangorestframework-csv rest-framework-generic-relations django-pdb munch py pytest python-dateutil svgwrite tabulate==0.7.5 Pillow
pip install matplotlib jupyter jupyterlab

# clean up
conda config --set auto_update_conda False
conda clean -y -a
apt-get clean -yq

%environment
export SINGULARITY_SHELL=/bin/shell.sh


