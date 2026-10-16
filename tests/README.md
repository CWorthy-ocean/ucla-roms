# ROMS pytest suite

## What the suite does

Each test compiles ROMS for one configuration, runs a few tens of time steps under MPI, joins the per-rank output with `ncjoin` and hashes every output variable with SHA-256. The hashes are compared with `results/results_<environ>.json`. Only `github_gnu` and `github_ifx` references exist, and only CI produces them. Laptop runs (`--environ=laptop`) have no reference file and are only useful relative to another local run.

Each variable's reference entry holds the hash and also the variable's min, max, mean and rms, rounded to 10 significant digits. The statistics are documentation: the comparison is on the hash alone. The reference files hold one line per variable, so the diff of a reference update shows each changed variable with its statistics.

## Why results are reproducible in CI

ROMS built at `-O0` is bitwise stable across GitHub runners. The inputs are not shipped: roms-tools generates them on the runner, and numpy, OpenBLAS and numba pick AVX-512 kernels when the VM exposes them (about one run in three), which changes the inputs at the last bit. The top-level `env:` block in `.github/workflows/containerized_ci.yml` pins all three to AVX2 paths. numpy's pin is the allow-list `NPY_ENABLE_CPU_FEATURES`, spelled with numpy's own dispatch-target names: numpy ignores a name it does not know with an ImportWarning that Python hides, and then dispatches as if the variable were unset (the group name `X86_V4` failed that way). `conftest.py` therefore fails the session when any target outside the allow-list is still enabled.

The `== input file hashes and environment ==` block printed at the end of every session is the sentinel: one SHA-256 over the stored variable values of each generated file (not the file bytes, which embed netCDF library versions). If the input hashes change without a Python-environment change (container rebuild, roms-tools or numpy bump), that is a bug in the pinning: investigate it before touching any reference. Each CI job also uploads the generated files themselves as the `inputs-<environ>` artifact (about 2 MB), so a sentinel change can be traced to the variables that moved: download two runs' artifacts and compare them with xarray. The block's `numpy_dispatch_enabled` list shows which SIMD targets numpy enabled; with the CI pinning, no `AVX512*` entry should appear on any runner.

## When a test fails on a hash mismatch

1. Read the failure message: it lists, per component (`physics`, `bgc`), each variable whose hash changed and the statistics that moved (`max 25.12345678 -> 25.12345679`). `statistics equal to 10 significant digits` means the change is either below about 1e-10 relative (roundoff) or moved values without changing min, max, mean or rms (a shifted or reordered field); the report of step 2 tells which. A mean close to zero, as in a symmetric field, can change in every digit from roundoff alone (e.g. `mean -5.7e-24 -> -1.7e-24`); judge it by the report's difference relative to the field's largest value.
2. Read the report on how far the outputs moved. The job reruns the suite on the base commit (the PR base, or the previous `main` commit on a push) and reports, for each changed variable, the number of points differing, the largest absolute difference and where it is, and that difference relative to the variable's largest base value. Read it on the run's summary page, or in the job log under "Report differences from the base commit". It runs only when the failure includes hash mismatches, and it takes a few extra minutes. The base run uses the base commit's own test code and input generation, so the report shows what the PR changed end to end.
3. Decide whether the change is intended by your PR. If not, it is a regression; fix the code.
4. If intended, take the run id from the failing CI run URL (`.../actions/runs/<run-id>`) and run, from anywhere in the checkout:
   ```
   python tests/update_references.py --run <run-id>
   ```
   It downloads the `results-github_gnu` and `results-github_ifx` artifacts with `gh`, merges them into `tests/results/results_<environ>.json` and prints what changed, separating `hash changed` from `statistics only` (hash equal, statistics moved, as after a numpy bump). Use `--dry-run` to look first.
5. Commit the two JSON files in the same PR and state the reason for the change in the PR body.

Both jobs upload their artifacts even when tests fail, so one failing run is enough to update both references. Artifacts are kept for 30 days.

## Adding or removing a test

A new test fails with "No reference results"; the same command adds it. For a removed test, run `update_references.py --run <run-id> --prune` on a CI run of the full suite: `--prune` deletes every reference key absent from the computed results, so it must never follow a `-k` subset run.

## Why not a bot

Most PRs come from a fork, where the workflow token is read-only and a bot cannot push regenerated references to the branch. The author runs the script locally instead.

## Running locally

```
export ROMS_ROOT=/path/to/ucla-roms
cd $ROMS_ROOT/tests
pytest -svvv -rfEx . --environ=github_gnu
```

Local hashes differ from the CI references because the compiler differs, so expect mismatches. To check a change locally, run with `--environ=laptop` before and after and compare the two printed result sets. `--tolerate-hash-mismatch` downgrades hash mismatches to XFAIL and exists for local exploration only; CI does not use it.

To see how far a change moved the output, as the CI report does, run the suite once in a checkout of the base commit and once in your branch, each with its own `ROMS_ROOT`, its own tools (`make -C Tools-Roms partit ncjoin COMPILER=gnu`, with that checkout's `Tools-Roms` first on `PATH`) and its own `--basetemp`. Run with `--tolerate-hash-mismatch`: local hashes never match the CI references, so every hashed test is recorded as a mismatch, and the mismatches file then lists them all.

```
cd $ROMS_ROOT/tests
pytest -svv -rfEx . --environ=github_gnu --tolerate-hash-mismatch --basetemp /tmp/pytest_base
```

Repeat in your branch with `--basetemp /tmp/pytest_head`, then compare from the branch checkout:

```
cd $ROMS_ROOT
python tests/compare_outputs.py --mismatches tests/results/mismatches_github_gnu.json --head-root /tmp/pytest_head --base-root /tmp/pytest_base
```

To compare any two joined output files, run `python tests/compare_outputs.py BASE.nc HEAD.nc`.
