#!/bin/bash

INPUT_DATA_DIR=$ROMS_ROOT/profiling/profiling_namelist/input_data_namelist

#NP_XI, NP_ETA iterables:
NP_X=(1 2 3 6) # 8 12 16 24 32)
NP_Y=(1 2 3 6) # 8 12 16 24 32)

# Create directory for runs
topdir=profiling_runs_NAMELIST_$(date +%Y%m%d_%H%M%S)
mkdir $topdir
cd $topdir && topdir=$(pwd) # get full path

for idx in "${!NP_X[@]}"; do
    cd $topdir

    npx=${NP_X[$idx]}
    npy=${NP_Y[$idx]}

    # Figure out how many nodes we need:
    ncpus=$(($npx*npy))
    if [ $((ncpus % 128)) == 0 ];then
	extra_node=0
    else
	extra_node=1
    fi
    nnodes=$(($ncpus/128 + $extra_node))

    echo "Run $idx: npx=$npx and npy=$npy ($ncpus CPUs, $nnodes nodes)"
    thisdir="${topdir}/run_${idx}_npx_${npx}_npy_${npy}"
    cp -rpv ../template_configuration_namelist/ ${thisdir}/
    cd ${thisdir} && thisdir=$(pwd) #Update var to use full path
    ln -s $ROMS_ROOT/Work/Makefile .
    mkdir $thisdir/output

    mkdir $thisdir/partitioned_inputs/
    cd $thisdir/partitioned_inputs/
    ln -s ${INPUT_DATA_DIR}/parent_grid.nc .
    ln -s ${INPUT_DATA_DIR}/parent_tides.nc .
    ln -s ${INPUT_DATA_DIR}/parent_surface_bulk_forcing.nc .
    ln -s ${INPUT_DATA_DIR}/parent_bgc_surface_forcing.nc .
    ln -s ${INPUT_DATA_DIR}/parent_phys_boundary_forcing.nc .
    ln -s ${INPUT_DATA_DIR}/parent_bgc_boundary_forcing.nc .
    ln -s ${INPUT_DATA_DIR}/parent_rst.20120201000000.nc .

    cd $thisdir/
    sed -i -e "s/TEMPLATE_NPX/${npx}/g" namelist.nml
    sed -i -e "s/TEMPLATE_NPY/${npy}/g" namelist.nml

    sed -i -e "s/TEMPLATE_NPX/${npx}/g" run_with_profiling_anvil.sh
    sed -i -e "s/TEMPLATE_NPY/${npy}/g" run_with_profiling_anvil.sh

    sed -i -e "s/TEMPLATE_NNODES/${nnodes}/g" run_with_profiling_anvil.sh
    sed -i -e "s/TEMPLATE_NCPUS/${ncpus}/g" run_with_profiling_anvil.sh

    sbatch run_with_profiling_anvil.sh
done
