ROMS settings
===============
Most ROMS settings are read at runtime from a namelist file, described :doc:`here <runtime_settings>`; start from the template at ``$ROMS_ROOT/src/namelist.nml``. Features of the model are enabled or disabled at compile time with C pre-processor macros, described :doc:`here <compile_time_settings>`. They are set in a ``cppdefs.opt`` file in your case directory, copied from ``$ROMS_ROOT/src/cppdefs.opt``.

.. toctree::
   :maxdepth: 1

   runtime_settings
   compile_time_settings
