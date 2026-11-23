Results change in every 3-D run with ``UV_VIS2`` and either a nonzero ``visc2`` or, with ``SPONGE``, a nonzero ``v_sponge`` (see Bug Fixes); runs without lateral viscosity are unchanged bit for bit.
- The change is largest in ``ubar``, ``vbar`` and ``zeta`` in the sponge layer next to open boundaries, which the sponge now damps.
