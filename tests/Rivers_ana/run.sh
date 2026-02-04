#!/bin/bash
#SBATCH --job-name=GoM_fwflux
#SBATCH --output=Test.out
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --account=ees250129
#SBATCH --export=ALL
#SBATCH --time=0:01:00

srun -n 6 roms benchmark.in
