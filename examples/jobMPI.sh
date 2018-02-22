#!/bin/bash
#SBATCH --partition=curso
#SBATCH --job-name=mpi_hello_world
#SBATCH -o out/mpi.%N.%j.out
#SBATCH -e error/slurm.%N.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --ntasks-per-node=6
#SBATCH --ntasks-per-core=1
#SBATCH --exclusive
#SBATCH --time=01:05 # time format (HH:MM)

cd $SLURM_SUBMIT_DIR
module load intel/tools-16.0.109
mpiexec.hydra -n 5 -ppn 5 ../ljmd-mpi.x < argon_108.inp
