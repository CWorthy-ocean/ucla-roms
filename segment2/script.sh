#!/bin/bash
#SBATCH --job-name=segment2
#SBATCH --output=/anvil/scratch/x-uheede/cstar/iceland1/tasks/segment2/logs/segment2.out
#SBATCH --partition=wholenode
#SBATCH --ntasks=256
#SBATCH --account=ees250129
#SBATCH --export=ALL
#SBATCH --mail-type=ALL
#SBATCH --time=48:00:00

set -e
ulimit -s unlimited
CSTAR_SLURM_MAX_WALLTIME='48:00:00' cstar blueprint run /anvil/scratch/x-uheede/cstar/_forge_bp_runs/Iceland1_256procs/blueprints/B_Iceland1_256procs.yaml --clobber --directives /anvil/scratch/x-uheede/cstar/iceland1/tasks/segment2/work/directives.yaml