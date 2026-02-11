Best Practices
==============

Model Output
------------

To avoid partial files, corrupt files and post-processing errors, we recommend adapting the following practice when setting output frequency of output files:

Boundary data frequency, 'bry_time' (1800 seconds default), should always be divisble by the timestep.

Restart file frequency should be equal to, or multiplier of history and diagnostics file frequency, including ext files when nesting. If your restart file frequency is daily, then history and ext file frequency should be daily as well.

When joining  partial files (including using extract_data_join), we recommend reading in only one file at a time, to avoid post-processing errors. I.e. avoid comands like: 'ncjoin *.nc' which is prone to errors.
