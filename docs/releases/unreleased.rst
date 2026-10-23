.. _unreleased:

Unreleased
----------

.. note::
    This release is currently in development

Breaking Changes
~~~~~~~~~~~~~~~~

- N/A

New Features
~~~~~~~~~~~~

- N/A

Bug Fixes
~~~~~~~~~

- N/A

Improvements
~~~~~~~~~~~~

- N/A

Miscellaneous
~~~~~~~~~~~~~

- CI: the ``test_gnu`` container runs with ``options: --shm-size=1g``. Transports and code paths are unchanged, so output hashes cannot move. ``test_ifx`` (OpenMPI) and the ``compile_*`` jobs are unchanged. The base-commit rerun step runs in the same container and gets the fix too. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
- CI: a new ``Print MPI and shared memory info`` step (``continue-on-error: true``) prints ``df -h /dev/shm``, ``ucx_info -v`` and ``mpichversion``. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
- CI: the ``/dev/shm`` sampler described above. It runs in the background of ``Run tests`` and is stopped by an ``EXIT`` trap, so pytest's exit status is unchanged. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
- CI: a new ``Show shared memory usage`` step (``if: always()``, ``continue-on-error: true``) prints the sampler log. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
- Tested locally in the image with stand-ins for pytest, run the way GitHub runs it (``bash -eo pipefail``). An unlinked 40 MB holder was caught with one ``ps`` snapshot, a failing exit status (3) came through unchanged, and the sampler stopped when the step ended. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
- Docs: ``environments/Containerfiles/README.md`` shows ``podman run --shm-size=1g`` and points to the workflow comment, which holds the full explanation. Podman's default is also 64 MB, so local users of the image can hit the crash too. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
- Open PRs pick up the fix only after merging ``main`` (or a new push), because rerunning an old run reuses its old workflow. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
