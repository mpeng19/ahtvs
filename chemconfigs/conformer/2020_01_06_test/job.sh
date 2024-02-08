#!/bin/bash
#SBATCH --partition={{jobspec.details.partitions|join(",")}}
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=4000

#source activate rdkit_2017_03_03
#source new-modules.sh
#module load legacy
#module load openbabel
#module load centos6/0.0.1-fasrc01
#module load legacy/0.0.1-fasrc01
#module load openbabel/2.3.2-fasrc01

# For running locally. If needing to run on cluster, write new script

INCHIKEY={{jobspec.inchikey}}
mv smiles.smi $INCHIKEY.smi

cwd=$(pwd)
#scfolder="/tmp/$(date +%Y%m%d%H%M%S%N)/"
#mkdir $scfolder
#cp *smi $scfolder
#cd $scfolder
confgenpath=$(python -c "import confgen, os; print(os.path.dirname(confgen.__file__))")
python $confgenpath/run_generator.py -g 8000 -e 0.30 -p 0.5 -E 3.5 -t 0.30 -n {{jobspec.details.conformer_count}} --fallback_to_align $INCHIKEY.smi

#cp * $cwd
#cd $cwd
rm -r $scfolder

