#!/bin/bash
# pulls prebuilt singularity containers https://github.com/verysure/miniconda3-rdkit-dftb
singularity pull -n rdkit-dftb.img shub://verysure/miniconda3-rdkit-dftb
# create and install confgen files to the overlay image
singularity image.create --size 2 confgen.img
sudo singularity exec --bind ./files:/files --overlay confgen.img rdkit-dftb.img /bin/bash /files/setup.sh
