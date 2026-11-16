Experimental, in development: the optional namelist group ``SL_ADV_SETTINGS`` (``sl_adv``, ``sl_c_lo``, ``sl_c_hi``, ``sl_cmax``; requires ``IMEX_IGW``) for swept-volume advection on faces with a large 2-D Courant number.
- So far ``sl_adv = .true.`` only computes the blend weights and reports them on an ``SL`` line; the solution does not change, and a namelist without the group runs as before.
