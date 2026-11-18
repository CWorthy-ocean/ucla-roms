Every namelist now needs the group ``IMEX_IGW_SETTINGS`` with the required switch ``imex_igw`` (no default), which states whether the run takes the IMEX internal-wave step; ROMS stops at startup, naming the key, when the group or the key is missing or the value does not match the build.
- A run without ``IMEX_IGW`` (most runs) adds ``&IMEX_IGW_SETTINGS imex_igw = .false. /`` to its namelist.
- A run built with ``IMEX_IGW`` adds ``imex_igw = .true.`` to its ``IMEX_IGW_SETTINGS`` group.
- Without ``IMEX_IGW`` the group's other keys are ignored (a startup note lists them), so the template's group runs in either build once ``imex_igw`` matches it.
- Namelists written by C-Star or Forge need the key too: C-Star needs a namelist tier for this release.
- In a build without ``IMEX_IGW`` an unknown or removed key in the group now stops the run, as in a build with it (the group used to be skipped).
