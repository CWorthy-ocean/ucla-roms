.. _unreleased:

Unreleased
----------

.. note::
    This release is currently in development

Breaking Changes
~~~~~~~~~~~~~~~~


- ``tests/results/results_<environ>.json`` changes shape again: each variable's entry holds its hash plus min, max, mean and rms; a branch that regenerated references before this lands must rerun ``tests/update_references.py`` after merging main. (`#392 <https://github.com/CWorthy-ocean/ucla-roms/pull/392>`_)

New Features
~~~~~~~~~~~~


- When pytest fails on hash mismatches, CI reruns the suite on the base commit and reports, per changed variable, the points that differ, the largest difference, where it is, and its size relative to the base field. (`#392 <https://github.com/CWorthy-ocean/ucla-roms/pull/392>`_)

  - The report appears on the run's summary page and in the job log; it runs only on hash-mismatch failures and skips workflow_dispatch runs.
  - ``tests/compare_outputs.py`` produces the same report for two local runs, or for any two output files.

- Reference entries carry each variable's min, max, mean and rms to 10 significant digits, and a mismatch message lists which statistics moved for each changed variable. (`#392 <https://github.com/CWorthy-ocean/ucla-roms/pull/392>`_)

  - Pass or fail is still decided by the hash alone; the statistics only describe the change.
  - Reference files hold one line per variable, so a reference update's diff shows each changed variable with its numbers.
  - ``update_references.py`` separates hash changes from statistics-only changes.

Bug Fixes
~~~~~~~~~

- N/A

Improvements
~~~~~~~~~~~~

- N/A

Miscellaneous
~~~~~~~~~~~~~

- ``tests/README.md`` explains how to read the statistics and the base-commit report, and how to compare two local runs. (`#392 <https://github.com/CWorthy-ocean/ucla-roms/pull/392>`_)
