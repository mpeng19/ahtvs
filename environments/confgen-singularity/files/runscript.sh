#!/bin/bash
. /opt/conda/etc/profile.d/conda.sh
conda activate chem
exec python /pythonlib/confgen/run_generator.py "$@"
