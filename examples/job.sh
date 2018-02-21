#!/bin/bash
#
# For more examples see:https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts
#
#SBATCH --partition=curso # request curso partition
#
#SBATCH --job-name=openmp_helloworld
#SBATCH -o out/counts.%N.%j.out # STDOUT %N is the node, %j is the job
#SBATCH -e error/slurm.%N.%j.err # STDERR
#
#SBATCH --exclusive           # ensure we get our nodes exclusively
#SBATCH --nodes=1             # run on one node
#SBATCH --ntasks=28           # we use 28 threads
#SBATCH --cpus-per-task=1     # ensure each thread is on a separate core
#SBATCH --ntasks-per-node=28
#
#SBATCH --time=30:00 # time format (HH:MM)
#
# change into directory where job was submitted
cd $SLURM_SUBMIT_DIR

# run OpenMP program using 28 threads
OMP_NUM_THREADS=2
.././ljmd-omp.x < argon_108.inp
#.././ljmd-serial.x < argon_108.inp
