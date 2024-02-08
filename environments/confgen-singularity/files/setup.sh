#!/bin/bash
mkdir -p /pythonlib
tar -xf /files/confgen.tar.gz -oC /pythonlib
echo 'export PYTHONPATH=/pythonlib' >> /environment
cat /files/runscript.sh > /.singularity.d/runscript 
chmod a+x /.singularity.d/runscript 
