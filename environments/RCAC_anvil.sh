module load gcc/11.2.0
module load openmpi/4.0.6
module load netcdf-fortran/4.5.3

export MPIHOME=${OPENMPI_HOME}/
export NETCDFHOME=${NETCDF}/
export PATH=${PATH}:${NETCDFHOME}/bin
export LIBRARY_PATH=${LIBRARY_PATH}:${NETCDFHOME}/lib

export ROMS_ROOT=$(cd ../ && pwd) #Set parent directory to this file as ROMS_ROOT
export PATH=$PATH:"$ROMS_ROOT/Tools-Roms"
