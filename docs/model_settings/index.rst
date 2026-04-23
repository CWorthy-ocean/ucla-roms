ROMS settings
===============
ROMS settings are now predominantly controlled at runtime via a namelist, which is described :doc:`here <runtime_settings>`. An example namelist is provided at ``$ROMS_ROOT/src/namelist.nml``. Certain features are enabled or disabled at compile time via C pre-processor macros, which are set in ``$ROMS_ROOT/src/cppdefs.opt``.

.. toctree::
   :maxdepth: 1

   runtime_settings
