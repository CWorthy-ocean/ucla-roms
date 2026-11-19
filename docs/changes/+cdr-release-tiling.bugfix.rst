Parameterized CDR releases spread over several cells (``cdr_hsc`` > 0) now give the same results to the last bit with any tiling (``NP_XI``, ``NP_ETA``). Their normalization was a sum of per-rank partial sums, whose rounding depended on the tiling.
- Results of runs with such releases can change at the rounding level.
