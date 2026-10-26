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
- ROMS now stops at startup when an averaging period is not a whole multiple of the MARBL timestep ROMS uses, which is ``marbl_timestep`` rounded down to a multiple of ``dt``. The check covers two settings: (`#400 <https://github.com/CWorthy-ocean/ucla-roms/pull/400>`_)

  - ``output_period_cdr``, when ``do_cdr_output`` and ``wrt_cdr_avg`` are on and ``cdr_monthly_averages`` is off
  - ``output_period_bgc_avg_dia``, when ``wrt_bgc_dia_avg`` is on

- With ``marbl_timestep > dt``, the averaged MARBL-held ``_cdr`` fields change: ``FG_CO2``, ``FG_ALT_CO2``, ``pCO2SURF``, ``pCO2SURF_ALT_CO2``, ``CO3``, ``CO3_ALT_CO2``, ``co3_sat_arag``, ``co3_sat_calc``, ``zsatarag``, ``zsatcalc``. Each MARBL update now counts toward the window of the MARBL step that computed it (see Bug Fixes). With ``marbl_timestep == dt`` the ``_cdr`` output is bit-for-bit identical. (`#400 <https://github.com/CWorthy-ocean/ucla-roms/pull/400>`_)
- ``_bgc_dia_avg`` records change in three ways: (`#400 <https://github.com/CWorthy-ocean/ucla-roms/pull/400>`_)

  - They are stamped at the end of each full ``output_period_bgc_avg_dia``, not one time step earlier.
  - They no longer include the diagnostics computed from the start-up state, so every averaged value changes slightly.
  - With ``marbl_timestep > dt``, MARBL diagnostics count toward the window of the MARBL step that computed them.

- With ``wrt_bgc_his`` off, ``_bgc_avg`` records likewise move to the end of each full period and no longer include the start-up state. (`#400 <https://github.com/CWorthy-ocean/ucla-roms/pull/400>`_)
- ``nrpf_bgc_avg_dia`` and ``nrpf_bgc_avg`` are now honored. Before, whenever ``wrt_bgc_his`` or ``wrt_bgc_avg`` was on, every averaged write opened a new file. ``_bgc_dia_avg`` and ``_bgc_avg`` now hold ``nrpf`` records per file, so tools that expect one record per file must read the record dimension. (`#400 <https://github.com/CWorthy-ocean/ucla-roms/pull/400>`_)

New Features
~~~~~~~~~~~~


- ``_uscl`` files carry ``avg_begin_time`` and ``avg_end_time`` (seconds, same units as ``ocean_time``), and a run whose length is not a multiple of ``output_period_uscl`` ends with a partial window instead of losing it. (`#398 <https://github.com/CWorthy-ocean/ucla-roms/pull/398>`_)
- ``marbl_config_file = ""`` runs MARBL with its default settings. (`#402 <https://github.com/CWorthy-ocean/ucla-roms/pull/402>`_)

Bug Fixes
~~~~~~~~~


- Upscaled boundary fluxes are no longer applied one window late in the parent, and the child's last window of export is no longer dropped, when the CDR forcing is built with ``Tools-Roms/uscl_to_cdr.py``. (`#398 <https://github.com/CWorthy-ocean/ucla-roms/pull/398>`_)
- ``cdr_output.F90``: after a restart, the first averaged ``_cdr`` record of the MARBL-held fields no longer differs from an uninterrupted run's. (`#400 <https://github.com/CWorthy-ocean/ucla-roms/pull/400>`_)

  - MARBL updates these fields only every ``marbl_timestep``. The averages counted each update for the steps after it, up to the next MARBL step. The restart file does not carry those values, so the restarted run used values recomputed from the restart state.
  - Each update now counts for the steps up to and including the MARBL step that produced it. MARBL applies its air-sea flux to DIC on that step, so each window's averaged flux matches its DIC change.
  - Tracers, pH and CDR sources were already bit-identical.

- ``bgc_io.F90``: ``_bgc_dia_avg`` records after a restart are now bit-identical to an uninterrupted run's, whether or not ``marbl_timestep`` equals ``dt``. (`#400 <https://github.com/CWorthy-ocean/ucla-roms/pull/400>`_)

  - The start-up call ``wrt_bgc(.true.)`` averaged the start-up state into the first window and advanced the window clock, because the guard was commented out. Every window was shifted one step and straddled each restart, where the start-up sample is recomputed from the restart state.
  - MARBL diagnostics use the same accumulation timing as ``_cdr``.

- ``_bgc_avg`` no longer averages the start-up state into its first window when ``wrt_bgc_his`` is off. (`#400 <https://github.com/CWorthy-ocean/ucla-roms/pull/400>`_)
- ``_bgc_dia_avg`` and ``_bgc_avg`` files hold ``nrpf_bgc_avg_dia`` / ``nrpf_bgc_avg`` records instead of one each. ``record_avg`` and ``record_dia_avg`` were reset on every call because the line clearing ``first_step_avg`` was commented out. (`#400 <https://github.com/CWorthy-ocean/ucla-roms/pull/400>`_)
- ``_uscl`` boundary exports (``ALK_add_*``, ``DIC_add_*`` and the CDR_LITE ``<tracer>_add_*``) now include horizontal tracer mixing across the open boundaries, not only advection, so a child's export closes its tracer budget. (`#399 <https://github.com/CWorthy-ocean/ucla-roms/pull/399>`_)
- The ``units`` attribute of the ``*_add_*`` variables now states a per-cell flux, tracer × m3/s (e.g. ``meq/m^3 * m3/s``), instead of ``meq/m^3s-1``. (`#399 <https://github.com/CWorthy-ocean/ucla-roms/pull/399>`_)
- File names too long for ROMS are no longer cut off silently; ROMS stops at startup naming the namelist key and the limit. (`#402 <https://github.com/CWorthy-ocean/ucla-roms/pull/402>`_)

  - A path of 256 characters or more in ``grdname``, ``inifile``, ``frcfiles``, ``cdr_file``, ``extract_file`` or ``marbl_config_file`` was cut off and then reported as a missing file.
  - An ``output_root_name`` or ``extract_root_name`` longer than 220 characters left output file names cut off, so files could overwrite each other.
  - CDR_LITE tracer and CDR gas-exchange output cut their file names at 99 characters, so with a root longer than about 74 characters each file overwrote the last.

- A missing MARBL settings file (``marbl_config_file``) now stops ROMS at startup; ROMS used to create it empty and run MARBL on its default settings without saying so. (`#402 <https://github.com/CWorthy-ocean/ucla-roms/pull/402>`_)
- MARBL output lists: the template's ``marbl_diagnostics_to_write`` default now writes the two diagnostics it names (it was read as one unrecognized name and wrote all 343), and full lists are no longer read past their end. (`#402 <https://github.com/CWorthy-ocean/ucla-roms/pull/402>`_)

  - Every unrecognized name in ``marbl_tracers_to_write`` is now reported, not only those before the first recognized one.


Improvements
~~~~~~~~~~~~


- ``Tools-Roms/uscl_to_cdr.py`` refuses ``_uscl`` files whose windows do not follow each other (e.g. a missing file) instead of converting them silently. (`#398 <https://github.com/CWorthy-ocean/ucla-roms/pull/398>`_)
- New precheck ``check_avg_ends_on_marbl_step`` (``precheck.F90``), next to ``check_output_divides_rst``. The error names the setting, the period and the MARBL timestep ROMS uses. (`#400 <https://github.com/CWorthy-ocean/ucla-roms/pull/400>`_)

  - The file-rollover precheck already requires the restart period to be a multiple of each output file's period. A restart therefore always falls on a MARBL step when averaging is on.

- PIO branches get the same changes. Each step does the same work as before, with one scalar check at startup. The cost should be negligible, but it is unverified at scale. (`#400 <https://github.com/CWorthy-ocean/ucla-roms/pull/400>`_)
- ``FORCING_FILES`` holds up to 10000 files (was 360); the list is sized to the files given, so the higher cap costs memory only while the namelist is read. (`#402 <https://github.com/CWorthy-ocean/ucla-roms/pull/402>`_)
- Namelist errors stop ROMS right after the namelist is read, before any value is used. (`#402 <https://github.com/CWorthy-ocean/ucla-roms/pull/402>`_)

Miscellaneous
~~~~~~~~~~~~~

- CI: the ``test_gnu`` container runs with ``options: --shm-size=1g``. Transports and code paths are unchanged, so output hashes cannot move. ``test_ifx`` (OpenMPI) and the ``compile_*`` jobs are unchanged. The base-commit rerun step runs in the same container and gets the fix too. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
- CI: a new ``Print MPI and shared memory info`` step (``continue-on-error: true``) prints ``df -h /dev/shm``, ``ucx_info -v`` and ``mpichversion``. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
- CI: the ``/dev/shm`` sampler described above. It runs in the background of ``Run tests`` and is stopped by an ``EXIT`` trap, so pytest's exit status is unchanged. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
- CI: a new ``Show shared memory usage`` step (``if: always()``, ``continue-on-error: true``) prints the sampler log. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
- Tested locally in the image with stand-ins for pytest, run the way GitHub runs it (``bash -eo pipefail``). An unlinked 40 MB holder was caught with one ``ps`` snapshot, a failing exit status (3) came through unchanged, and the sampler stopped when the step ended. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
- Docs: ``environments/Containerfiles/README.md`` shows ``podman run --shm-size=1g`` and points to the workflow comment, which holds the full explanation. Podman's default is also 64 MB, so local users of the image can hit the crash too. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
- Open PRs pick up the fix only after merging ``main`` (or a new push), because rerunning an old run reuses its old workflow. (`#401 <https://github.com/CWorthy-ocean/ucla-roms/pull/401>`_)
- ``docs/model_settings/runtime_settings.md``: ``output_period_cdr`` and ``output_period_bgc_avg_dia`` state the multiple-of-MARBL-timestep requirement, and ``output_period_bgc_avg`` says what each record averages. (`#400 <https://github.com/CWorthy-ocean/ucla-roms/pull/400>`_)
- ``src/Make.depend`` is regenerated, since ``precheck`` now uses ``scalars`` and ``marbl_driver``. The other line changes there are drift that was already on ``main``. (`#400 <https://github.com/CWorthy-ocean/ucla-roms/pull/400>`_)
- New test: the interior inventory change of an ALK perturbation in a SPONGE + UPSCALING run matches minus the time-integrated ``_uscl`` export to round-off. (`#399 <https://github.com/CWorthy-ocean/ucla-roms/pull/399>`_)
- Corrected the sponge-profile comment in ``set_nudgcof`` and the ``UPSCALING`` description in ``cppdefs.opt`` and the compile-time settings docs. (`#399 <https://github.com/CWorthy-ocean/ucla-roms/pull/399>`_)
- Every file-name variable and its broadcast use ``max_name_size`` (now in ``param``), so raising it is a one-line change. (`#402 <https://github.com/CWorthy-ocean/ucla-roms/pull/402>`_)
- The test namelist helper asks for MARBL's default settings and every diagnostic, which is what the tests ran with before; new ``TestFileNames`` tests cover these limits. (`#402 <https://github.com/CWorthy-ocean/ucla-roms/pull/402>`_)
