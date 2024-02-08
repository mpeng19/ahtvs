#!/bin/bash
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -t 3-0:00 
#SBATCH -p {{jobspec.details.partitions|join(",")}}
#SBATCH --mem-per-cpu=5250

# Settings
INPUT_FILE="gaussian_test.inp"
OUTPUT_FILE="gaussian_test.log"
CUR_DIR=$(pwd)
SCRATCH=/scratch/
RESULTS=$(pwd)
G09_DIR=/n/sw/g09_D.01/g09/

module load centos6/0.0.1-fasrc01
module load gaussian/09_D.01-fasrc01



# Add blank line to end of file just in case
echo "" >> gaussian_test.inp

time_run(){

{ time -p g09 $INPUT_FILE > $OUTPUT_FILE & pid=$!; } 2>> $OUTPUT_FILE

}

if [ -d $SCRATCH ] && [[ $1 != noscratch ]]; then
    #setup scratch dirs
    WORK_DIR=$SCRATCH/$USER/$SLURM_JOB_ID
    echo "Setting up directories"
    echo "  at ${WORK_DIR}"
    mkdir -p $WORK_DIR
    cp -R $CUR_DIR/* $WORK_DIR
    cd $WORK_DIR

    clean_up()
    {
    # move results of calc
    touch job_manager-complete
#    echo "Copying results back to ${RESULTS}"
#    rsync -a $WORK_DIR/* $RESULTS/ --exclude=slurm-*
    cp $WORK_DIR/job_manager-complete $RESULTS/
    cp $WORK_DIR/*.log $RESULTS/
    cp $WORK_DIR/*.chk $RESULTS/
    echo 'Cleaning up'
    rm -R $WORK_DIR/*
    cd $CUR_DIR
    }

    trap 'kill -9 $pid; clean_up; exit' SIGTERM SIGINT

    time_run

    wait
    clean_up
else
    time_run

    wait
fi




