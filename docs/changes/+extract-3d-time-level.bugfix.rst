Extracted data, including the boundary files written for nested child domains, now hold ``u``, ``v``, ``temp``, ``salt`` and the other tracers at the time of their record, like ``zeta``, ``ubar`` and ``vbar``; they were one parent time step older.
- Extracted 3D fields change accordingly: each record now holds what the previous code wrote in the next record.
