Singularity Container
=====================

Singularity containter with ubuntu (LTS), miniconda3, rdkit, django and jupyter.
The `latest` (default) builds with latest os, conda and python packages.
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
* [Django](https://www.djangoproject.com/)
* [scipy](https://www.scipy.org/) stack
* [Jupyter](https://jupyter.org/) 
