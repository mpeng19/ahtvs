#!/bin/bash
#SBATCH --partition={{jobspec.details.partitions|join(",")}}
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 75:00:00
#SBATCH --mem-per-cpu=4000

source new-modules.sh
module load legacy
module load openbabel
module load hpc/mongoengine-10.17.2013
module load hpc/python-2.7.3
module load centos6/rdkit-2013.09.1_with_inchi

source activate a2g2

INCHIKEY={{jobspec.inchikey}}
mv smiles.smi $INCHIKEY.smi

cwd=$(pwd)
scfolder="/tmp/$(date +%Y%m%d%H%M%S%N)/"
mkdir $scfolder
cp *smi $scfolder
cd $scfolder
confgenpath=$(python -c "import confgen, os; print(os.path.dirname(confgen.__file__))")
python $confgenpath/run_generator.py -g 8000 -e 0.30 -p 0.5 -E 3.5 -t 0.30 -n {{jobspec.details.conformer_count}} --fallback_to_align $INCHIKEY.smi

cp * $cwd
cd $cwd
rm -r $scfolder

