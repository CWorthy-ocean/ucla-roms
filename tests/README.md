# ROMS pytest suite

## What the suite does

Each test compiles ROMS for one configuration, runs a few tens of time steps under MPI, joins the per-rank output with `ncjoin` and hashes every output variable with SHA-256. The hashes are compared with `results/results_<environ>.json`. Only `github_gnu` and `github_ifx` references exist, and only CI produces them. Laptop runs (`--environ=laptop`) have no reference file and are only useful relative to another local run.

## Why results are reproducible in CI

ROMS built at `-O0` is bitwise stable across GitHub runners. The inputs are not shipped: roms-tools generates them on the runner, and numpy, OpenBLAS and numba pick AVX-512 kernels when the VM exposes them (about one run in three), which changes the inputs at the last bit. The top-level `env:` block in `.github/workflows/containerized_ci.yml` pins all three to AVX2-generic paths.

The `== input file hashes and environment ==` block printed at the end of every session is the sentinel: one SHA-256 over the stored variable values of each generated file (not the file bytes, which embed netCDF library versions). If the input hashes change without a Python-environment change (container rebuild, roms-tools or numpy bump), that is a bug in the pinning: investigate it before touching any reference.

## When a test fails on a hash mismatch

1. Read the failure message: it lists which variables changed, per component (`physics`, `bgc`).
2. Decide whether the change is intended by your PR. If not, it is a regression; fix the code.
3. If intended, take the run id from the failing CI run URL (`.../actions/runs/<run-id>`) and run, from anywhere in the checkout:
   ```
   python tests/update_references.py --run <run-id>
   ```
   It downloads the `results-github_gnu` and `results-github_ifx` artifacts with `gh`, merges them into `tests/results/results_<environ>.json` and prints what changed. Use `--dry-run` to look first.
4. Commit the two JSON files in the same PR and state the reason for the change in the PR body.

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
