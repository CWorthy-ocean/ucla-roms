Nested simulation with carbon dioxide removal
=============================================

This tutorial describes a ROMS set-up for conducting a carbon dioxide removal (CDR) intervention simulation via ocean alkalinity enhancement (OAE).
It consists of three ROMS simulations:

1. a :doc:`lower-resolution "parent" simulation <1_parent_simulation>` of the northeast Atlantic, which generates initial and boundary conditions for...
2. a :doc:`higher-resolution "child" simulation <2_child_simulation_with_cdr>` of the northwest coast of Iceland, where an alkalinity perbuation is made and eventually fluxes out into...
3. a :doc:`second run of the parent simulation <3_parent_simulation_re-run_with_upscaling>`, which tracks CDR-relevant tracers that have left the child domain in (2).

.. toctree::
   :maxdepth: 1

   1. Parent simulation <1_parent_simulation>
   2. Child simulation with CDR <2_child_simulation_with_cdr>
   3. Parent simulation re-run with upscaled CDR <3_parent_simulation_re-run_with_upscaling>
