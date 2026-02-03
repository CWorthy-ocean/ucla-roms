#!/bin/bash
# Script to download input data used by Example cases
if command -v curl &> /dev/null;then
    DOWNLOAD_COMMAND="curl -O -L"
elif command -v wget &> /dev/null;then
    DOWNLOAD_COMMAND=wget
fi
echo ${DOWNLOAD_COMMAND}
# URLs to download
URL_PREFIX="https://github.com/CWorthy-ocean/ucla_roms_example_input_data/raw/main"
files=(
    "example_input_surface_forcing.nc"
    "example_input_pipe_forcing.nc"
    "example_input_bgc_surface_forcing_clim.nc"
    "example_input_tides.nc"
    "example_input_river_forcing.nc"
    "example_input_boundary_forcing.nc"
    "example_input_bgc_boundary_forcing_clim.nc"
    "example_input_surface_flux_forcing.nc"
    "example_input_grid.nc"
    "example_input_bgc_initial_conditions.nc"
    "cdr_forcing_3d.nc"
)

for fname in "${files[@]}";do
    echo "#######################################################"
    echo "FETCHING FILE ${fname}"
    echo "#######################################################"
    ${DOWNLOAD_COMMAND} "${URL_PREFIX}/${fname}"
    partit 3 2 "${fname}"
    #echo "${URL_PREFIX}/${fname}"
done

# Unpartitioned CDR forcing:
${DOWNLOAD_COMMAND} "${URL_PREFIX}/cdr_forcing_dp.nc"
${DOWNLOAD_COMMAND} "${URL_PREFIX}/cdr_forcing_parm.nc"
