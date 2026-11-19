With ``pio_direct_read``, initial-condition files written by roms-tools (whose ``ocean_time`` dimension is fixed, not unlimited) are now read directly too; before, all their variables went through PIO.
- Rank 0 now says, once per file, why any variable was read through PIO: how many, the first names and the reason.
