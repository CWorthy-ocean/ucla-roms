#!/usr/bin/env python3
"""
cppcheck.py - Scan cppdefs.opt for CPP macro names and generate
check_switches1.F, a Fortran subroutine that records which switches
were active at compile time.

This is a Python rewrite of the original cppcheck.F. The "main
configuration" guard logic has been removed; any switch that appears
in cppdefs.opt is cataloged unconditionally.
"""

import re
import sys
from pathlib import Path

SOURCE = "cppdefs.opt"
OUTPUT = "check_switches1.F"

# CPP directives that introduce macro names we want to catalog.
# After one of these, the next identifier is a macro name.
MACRO_INTRODUCERS = {"define", "undef", "ifdef", "ifndef", "if",
                     "elif", "defined"}

# Other directives — listed so we don't mistake them for macro names.
OTHER_DIRECTIVES = {"else", "endif", "include", "pragma", "error",
                    "warning", "line"}


def extract_switches(path: Path) -> list[str]:
    """Return CPP macro names found in `path`, in order, deduplicated."""
    text = re.sub(r"/\*.*?\*/", " ", path.read_text(), flags=re.DOTALL)
    seen: set[str] = set()
    switches: list[str] = []

    for raw_line in text.splitlines():
        if not raw_line.lstrip().startswith("#"):
            continue
        tokens = re.findall(r"[A-Za-z_][A-Za-z0-9_]*", raw_line)

        i = 0
        while i < len(tokens):
            if tokens[i] in MACRO_INTRODUCERS and i + 1 < len(tokens):
                name = tokens[i + 1]
                if (name not in MACRO_INTRODUCERS
                        and name not in OTHER_DIRECTIVES
                        and name not in seen):
                    seen.add(name)
                    switches.append(name)
                i += 2
            else:
                i += 1

    return switches


HEADER = f'''\
#include "cppdefs.opt"

      subroutine check_switches1(ierr)

! WARNING: machine-generated file. Do not edit by hand.
! Regenerate with cppcheck.py whenever a new CPP-switch is
! introduced into "{SOURCE}". No action is needed if changes
! are limited to activating or deactivating previously-known
! switches.
! Total number of CPP-switches: {{n}}

      use param
      use scalars
      use strings
      use error_handling_mod, only: error_log
      implicit none
      integer ierr, is, ie

      do is=1,max_opt_size
        cpps(is:is)=' '
      enddo

      is=1 ; ie=is+{len(SOURCE) + 1}
      if (ie >= max_opt_size) call insufficient_size
      cpps(is:ie)='<{SOURCE}>'
      cpps(ie+1:ie+1)=' '
'''

SWITCH_BLOCK = '''\
#ifdef {name}
      is=ie+2 ; ie=is+{end}
      if (ie >= max_opt_size) call insufficient_size
      cpps(is:ie)='{name}'
      cpps(ie+1:ie+1)=' '
#endif
'''

FOOTER = '''\
      call error_log%abort_check()

      return

      contains

      subroutine insufficient_size
      call error_log%raise_global(context="check_switches1",
     &  info="Insufficient size of string 'cpps' in 'strings.F'."//
     &       " Increase its size and recompile")
      end subroutine insufficient_size

      end
'''

def generate_fortran(switches: list[str]) -> str:
    blocks = [SWITCH_BLOCK.format(name=name, end=len(name) - 1)
              for name in switches]
    return HEADER.format(n=len(switches)) + "".join(blocks) + FOOTER

def main() -> int:
    src = Path(SOURCE)
    if not src.is_file():
        print(f'### ERROR: cannot read {src}', file=sys.stderr)
        return 1

    print(f'cppcheck.py: scanning {src}')
    switches = extract_switches(src)
    print(f'  found {len(switches)} CPP-switches')
    for name in switches:
        print(f'    {name}')

    Path(OUTPUT).write_text(generate_fortran(switches))
    print(f'  wrote {OUTPUT}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
