.. _unreleased:

Unreleased
----------

.. note::
    This release is currently in development

Breaking Changes
~~~~~~~~~~~~~~~~


- ``_uscl`` files no longer have a start-up record at the run start: a run of ``n`` whole windows writes ``n`` records, not ``n + 1``, and the first file is created at the run start but written at the end of the first window. (`#398 <https://github.com/CWorthy-ocean/ucla-roms/pull/398>`_)
- ``_uscl`` window averages now cover the steps of the window ending at ``ocean_time`` exactly; before, each average was shifted one step late and the last window counted its last step twice, so ``_add_`` and ``h_`` values change slightly. (`#398 <https://github.com/CWorthy-ocean/ucla-roms/pull/398>`_)
- ``Tools-Roms/uscl_to_cdr.py`` stamps each CDR record at the start of its window and appends a zero-flux record at the end of the last window, so its output has one more record than the ``_uscl`` windows and different ``cdr_time`` values. (`#398 <https://github.com/CWorthy-ocean/ucla-roms/pull/398>`_)

New Features
~~~~~~~~~~~~


- ``_uscl`` files carry ``avg_begin_time`` and ``avg_end_time`` (seconds, same units as ``ocean_time``), and a run whose length is not a multiple of ``output_period_uscl`` ends with a partial window instead of losing it. (`#398 <https://github.com/CWorthy-ocean/ucla-roms/pull/398>`_)

Bug Fixes
~~~~~~~~~


- Upscaled boundary fluxes are no longer applied one window late in the parent, and the child's last window of export is no longer dropped, when the CDR forcing is built with ``Tools-Roms/uscl_to_cdr.py``. (`#398 <https://github.com/CWorthy-ocean/ucla-roms/pull/398>`_)

Improvements
~~~~~~~~~~~~


- ``Tools-Roms/uscl_to_cdr.py`` refuses ``_uscl`` files whose windows do not follow each other (e.g. a missing file) instead of converting them silently. (`#398 <https://github.com/CWorthy-ocean/ucla-roms/pull/398>`_)

Miscellaneous
~~~~~~~~~~~~~

- CI: the ``test_gnu`` container runs with ``options: --shm-size=1g``. Transports and code paths are unchanged, so output hashes cannot move. ``test_ifx`` (OpenMPI) and the ``compile_*`` jobs are unchanged. The base-commit rerun step runs in the same container and gets the fix too. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
- CI: a new ``Print MPI and shared memory info`` step (``continue-on-error: true``) prints ``df -h /dev/shm``, ``ucx_info -v`` and ``mpichversion``. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
- CI: the ``/dev/shm`` sampler described above. It runs in the background of ``Run tests`` and is stopped by an ``EXIT`` trap, so pytest's exit status is unchanged. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
- CI: a new ``Show shared memory usage`` step (``if: always()``, ``continue-on-error: true``) prints the sampler log. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
- Tested locally in the image with stand-ins for pytest, run the way GitHub runs it (``bash -eo pipefail``). An unlinked 40 MB holder was caught with one ``ps`` snapshot, a failing exit status (3) came through unchanged, and the sampler stopped when the step ended. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
- Docs: ``environments/Containerfiles/README.md`` shows ``podman run --shm-size=1g`` and points to the workflow comment, which holds the full explanation. Podman's default is also 64 MB, so local users of the image can hit the crash too. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
- Open PRs pick up the fix only after merging ``main`` (or a new push), because rerunning an old run reuses its old workflow. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
