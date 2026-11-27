Every namelist now needs the group ``IMEX_IGW_SETTINGS`` with the required switch ``imex_igw`` (no default), which states whether the run takes the IMEX internal-wave step; ROMS stops at startup, naming the key, when the group or the key is missing.
- A run without IMEX (most runs) adds ``&IMEX_IGW_SETTINGS imex_igw = .false. /`` to its namelist.
- With ``imex_igw = .false.`` the group's other keys are ignored (a startup note lists them), so the template's group runs either way.
- Namelists written by C-Star or Forge need the key too: C-Star needs a namelist tier for this release.
- An unknown or removed key in the group stops the run, also with ``imex_igw = .false.``.
- The IMEX step was first selected with a cpp key ``IMEX_IGW`` that ``imex_igw`` had to match; that key is gone, and a ``cppdefs.opt`` that defines it stops the build.
