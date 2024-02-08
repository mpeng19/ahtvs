#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -t 1400
#SBATCH --partition=aspuru-guzik,serial_requeue
#SBATCH --mem-per-cpu=2000

# Settings
INPUT_FILE="bp86_6-31gs_opt.inp"
OUTPUT_FILE="bp86_6-31gs_opt.out"

source /etc/profile.d/modules.sh

clean_up()
{
echo 'Cleaning up'
scratch_folder_number=`grep scratch slurm*.out | grep -oE '[^ ]+$' | grep -oE '[0-9]{1,7}'`
echo 'Folder number is' qchem"$scratch_folder_number"
user=`whoami`
if  [ $scratch_folder_number ]
    then
       if [ -e  /scratch/qchem"$scratch_folder_number"'.0'  ]
            then
            find /scratch/ -maxdepth 1  -type d -user $user  -name qchem"$scratch_folder_number""*"  -exec rm -rf {} \;
            echo 'removed /scratch/qchem'"$scratch_folder"' folders'
            if [ -e  /scratch/qchem"$scratch_folder"'.0'  ]
                then
                echo 'Folder not removed - something went wrong'
            fi
       fi
fi
}

trap 'kill -9 $pid; clean_up; exit' SIGTERM SIGINT

module load centos6/qchem-4.2
qchem -nt 8 $INPUT_FILE $OUTPUT_FILE & pid=$!
wait
clean_up
