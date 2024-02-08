[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/980)


Singularity Container
=====================

Singularity containter with ubuntu (LTS), miniconda3, rdkit and dftb.
The `latest` (default) builds with latest os, conda and python packages.
The `stable` uses ubuntu 18.04, conda3 4.5.2 and other freeze python
packages.
It also contains custom shell entry for sourcing `conda.sh`.


Building Manually
-----------------

To build the image, run `sudo singularity build <name.img> Singularity`.
See [Singularity](https://singularity.lbl.gov/) 
for more info. 


Acknowledge
-----------

The recipes build from many open source projects, including
* [Singularity](https://singularity.lbl.gov/)
* [Ubuntu](https://www.ubuntu.com/)
* [Conda](https://conda.io/)
* [RDkit](http://www.rdkit.org/)
* [DFTB](http://www.dftb.org/) and [DFTB+](http://www.dftbplus.org/)
* [scipy](https://www.scipy.org/) stack

