Start-up reads the time axes of classic-format (CDF-1, CDF-2, CDF-5) forcing files in 8-KiB blocks instead of one 4-MiB block per time value, so a run no longer reads gigabytes to find its first forcing records; results are unchanged.
- On a Pacific MARBL run on 2 nodes from a cold cache, the forcing search before the first record drops from 29-36 s to 7-8 s.
