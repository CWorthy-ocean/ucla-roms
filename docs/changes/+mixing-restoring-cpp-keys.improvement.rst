Lateral viscosity and tracer diffusion are skipped when they are zero everywhere, which saves their cost in runs with ``visc2`` and ``tnu2`` of 0 and no boundary sponge.
- ``tnu2`` (``TRACER_DIFF2``) and ``akt_bak`` (``VERTICAL_MIXING_SETTINGS``) values beyond the number of tracers are ignored, instead of stopping the run.
