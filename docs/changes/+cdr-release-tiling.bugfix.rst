The normalization of parameterized CDR releases spread over several cells (``cdr_hsc`` > 0) no longer depends on the tiling (``NP_XI``, ``NP_ETA``): it was a sum of per-rank partial sums, whose rounding depended on the tiling.
- With ``BUILD_MODE=strict`` or ``BUILD_MODE=test``, runs with such releases now give the same results to the last bit with any tiling. The default optimized build has other tiling-dependent rounding, with or without CDR releases, but its release weights are the same in every tiling.
- Results of runs with such releases can change at the rounding level.
