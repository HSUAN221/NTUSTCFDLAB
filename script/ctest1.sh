#!/bin/bash
#PBS -P MST108362
#PBS -N atest_0012
#PBS -l select=1:ncpus=20:mpiprocs=5:ompthreads=4
#PBS -l walltime=00:30:00
#PBS -q ctest
#PBS -j oe
#PBS -m be 








module purge
module load intel/2018_u1
module list
cd $PBS_O_WORKDIR


cat $PBS_NODEFILE
echo $PBS_O_WORKDIR
date

mpirun -PSM2 ./sol0


#qstat -xs
 
