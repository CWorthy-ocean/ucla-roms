#!/bin/bash

#SBATCH -A ees250129  # Allocation name
#SBATCH --nodes=TEMPLATE_NNODES        # Total # of nodes 
#SBATCH --ntasks=TEMPLATE_NCPUS     # Total # of MPI tasks
#SBATCH --time=00:05:00   # Total run time limit (hh:mm:ss)
#SBATCH -p wholenode 
#SBATCH -J ROMS_profile_TEMPLATE_NCPUS     # Job name
#SBATCH -o roms.out     # Name of stdout output file
#SBATCH --mail-user=sam.maticka@cworthy.org
#SBATCH--mail-type=all   # Send email to above address at begin and end of job

set -e

# Partition inputs
cd partitioned_inputs
start_time=$(date +%s)
for F in *.nc;do
    partit TEMPLATE_NPX TEMPLATE_NPY ${F};
done
end_time=$(date +%s)
cd ../
partitioning_time=$((end_time - start_time))

# Compile
source /home/x-smaticka/.ROMS
start_time=$(date +%s)
make PROFILER=scorep
end_time=$(date +%s)
compile_time=$((end_time - start_time))

# Run
start_time=$(date +%s)
srun -n TEMPLATE_NCPUS ./roms ./namelist.nml
end_time=$(date +%s)
run_time=$((end_time - start_time))

# Joint outputs
ncpus=TEMPLATE_NCPUS
nchar=$(echo ${ncpus} | wc -c)
QM=$(printf "%0.s?" $(seq 1 $((nchar-1)))) #Sequence of ?s length of ncpus
ZR=$(printf "%0.s0" $(seq 1 $((nchar-1)))) #Sequence of 0s length of ncpus

cd output
start_time=$(date +%s)
for F in *.${ZR}.nc; do
    if [[ "${F}" == *"_ext."* ]]; then
        continue
    fi    
    ncjoin ${F/.${ZR}.nc}.${QM}.nc
done
end_time=$(date +%s)
join_time=$((end_time - start_time))

echo "SIMULATION COMPLETE"
echo "PARTITIONING TIME: ${partitioning_time}s"
echo "COMPILE TIME: ${compile_time}s"
echo "RUN TIME: ${run_time}s"
echo "JOIN TIME: ${join_time}s"
