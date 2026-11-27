``imex_corr_adv_vimp = 1`` in ``IMEX_IGW_SETTINGS`` gives an IMEX cell whose vertical Courant number exceeds 1 the corrector's advective buoyancy change (``imex_corr_adv = 1``) from the tracer step's own estimate, explicit with ``We`` and then the implicit upwind step with ``Wi``, in place of no change.
- The default 0 keeps no change in those cells, and the results are unchanged bit for bit.
- Swept-volume vertical advection needs 1.
