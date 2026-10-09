.. _unreleased:

Unreleased
----------

.. note::
    This release is currently in development

Breaking Changes
~~~~~~~~~~~~~~~~


- ``tests/results/results_<environ>.json`` changes shape from one digest per output file to one digest per variable; existing reference files are regenerated in this PR, and any out-of-tree reference file must be regenerated with ``tests/update_references.py``. (`#384 <https://github.com/CWorthy-ocean/ucla-roms/pull/384>`_)

New Features
~~~~~~~~~~~~


- ``python tests/update_references.py --run <run-id>`` downloads a CI run's results artifacts and merges them into the reference files, printing which tests and variables changed. (`#384 <https://github.com/CWorthy-ocean/ucla-roms/pull/384>`_)

  - ``--from-file`` merges a local computed file, ``--prune`` drops references for removed tests after a full-suite run, ``--dry-run`` reports without writing.

- A hash mismatch now names the variables that changed, were added or were removed, and a missing reference is reported separately from a mismatch; both messages print the update command. (`#384 <https://github.com/CWorthy-ocean/ucla-roms/pull/384>`_)

Bug Fixes
~~~~~~~~~

- N/A

Improvements
~~~~~~~~~~~~


- Each pytest session records the SHA-256 of every generated input file's variable values plus the roms-tools, numpy, scipy, numba and xarray versions, the runner's AVX-512 exposure and the SIMD dispatch targets numpy enabled, so an input change can be told apart from a model change. (`#384 <https://github.com/CWorthy-ocean/ucla-roms/pull/384>`_)
- The session fails at start if ``NPY_ENABLE_CPU_FEATURES`` is set but numpy still enables a target outside it, so an ineffective dispatch pin is reported instead of silently producing runner-dependent inputs. (`#384 <https://github.com/CWorthy-ocean/ucla-roms/pull/384>`_)
- ``--inputs-out DIR`` copies the generated (unpartitioned) input files to a directory; CI uploads them as the ``inputs-github_gnu`` / ``inputs-github_ifx`` artifacts so a sentinel change can be traced to the variables that moved. (`#384 <https://github.com/CWorthy-ocean/ucla-roms/pull/384>`_)
- ``--results-out`` writes the computed results JSON at session end (default ``tests/results/computed_<environ>.json``), replacing copy-from-log reference updates. (`#384 <https://github.com/CWorthy-ocean/ucla-roms/pull/384>`_)

Miscellaneous
~~~~~~~~~~~~~

- CI pins numpy, OpenBLAS and numba to AVX2 code paths so roms-tools generates identical inputs on every runner, and no longer passes ``--tolerate-hash-mismatch``, so a reference mismatch fails the job again. (`#384 <https://github.com/CWorthy-ocean/ucla-roms/pull/384>`_)
- Both pytest jobs upload ``results-github_gnu`` / ``results-github_ifx`` artifacts with the computed results and input sentinel, retained 30 days, also on failure. (`#384 <https://github.com/CWorthy-ocean/ucla-roms/pull/384>`_)
- ``tests/README.md`` documents how the suite works, the reference-update loop, adding and removing tests, and local runs. (`#384 <https://github.com/CWorthy-ocean/ucla-roms/pull/384>`_)
