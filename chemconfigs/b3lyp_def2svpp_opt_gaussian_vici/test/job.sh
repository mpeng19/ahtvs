#!/bin/bash
#
#PBS -l walltime=48:00:00,mem=4gb,nodes=1:ppn=1

# clean up at the end of the job
function cleanup {
 
 echo -n END DATE: 
 date 
 pwd

 \cp --backup=t -p $TMPDIR/gaussian_test.log $PBS_O_WORKDIR
   
 echo contents of $TMPDIR dir
 ls -al
 
 }
trap cleanup 0 1 2 3 9 15

module purge

 export g16root=/share/apps/Gaussian/g16_b01
 . $g16root/g16/bsd/g16.profile

cd $TMPDIR

env

echo -n START DATE: 
date 


g16  <  $PBS_O_WORKDIR/gaussian_test.inp  >  gaussian_test.log

exit
