Best Practices
==============

Model Output
------------

To avoid partial files, corrupt files and post-processing errors, we recommend adapting the following practice when setting output frequency of output files:

Output frequency should always be divisble by the timestep. Note, if you are nesting and have set ´do_extract = .true.´ in ´extract_data.opt´, your most frequent output is likely to be 1800 seconds, which is the default output frequency of boundary data (.ext). So your timestep should be divisible by 1800. 

Restart output frequency should be equal to, or multiplier of your history and diagnostics output file frequency, including ´.ext´ files when nesting to avoid partial files. If your restart file frequency is daily, then you should have maximum one history/diagnostics/cstar file written out per 24 hours. So, for example for history files, set your restart output period so that ´output_period_rst´=´output_period_his´*´nrpf_his´ / a, where a is an integer of at least one.

When joining  partial files (including using extract_data_join), we recommend reading in only one file at a time, to avoid post-processing errors. I.e. avoid comands like: 'ncjoin *.nc' which is prone to errors.
