"""Helpers for asserting test output hashes match the reference JSON.

Every call to :func:`assert_output_matches_reference` records its computed
values into the module-level :data:`computed_results` dict, regardless of
whether the assertion passed or failed. ``conftest.py``'s
``pytest_sessionfinish`` hook then writes the collected dict to a JSON file at
the end of the run, which ``tests/update_references.py`` turns into
``results/results_<environ>.json`` to update references.

Each variable's entry holds its SHA-256 and some statistics of its values. Only
the hash decides whether a test passes; the statistics say by how much a
mismatching variable moved. Every mismatch is also recorded in
:data:`mismatched_outputs`, which names the output files that disagree.
"""
from typing import TYPE_CHECKING

import pytest

from .update_references import predates_statistics

# _helpers reads ROMS_ROOT from the environment when it is imported, so it is
# imported for the type checker only.
if TYPE_CHECKING:
    from ._helpers import OutputSummary

# Session-wide collector, populated by `assert_output_matches_reference` and
# read by the `pytest_sessionfinish` hook in conftest.py. Only tests that
# actually ran will appear here, so running a subset of the suite produces a
# partial dict rather than entries with missing values.
computed_results: dict[str, dict] = {}

# Set from the `--environ` command-line flag by `pytest_configure` in
# conftest.py. It lives here, like `tolerate_hash_mismatch`, because the
# assertion helper has no access to the pytest config; it is used to name the
# reference file in failure messages.
environ: str = "laptop"

# SHA-256 of every generated input file plus the Python environment and CPU
# feature that produced them. Filled by the `input_dir` fixture in conftest.py
# and written out by `pytest_sessionfinish`. The inputs are generated on the
# runner and have varied with CPU features, so this record is what tells an
# input change apart from a ROMS change when output hashes move.
input_file_hashes: dict = {}

# Set from the `--tolerate-hash-mismatch` command-line flag by
# `pytest_configure` in conftest.py. When true, a hash mismatch is reported as
# an expected failure (XFAIL) instead of a failure, so the run stays green.
# Every other way a test can fail -- a compile error, a ROMS crash, a missing
# output file -- is untouched and still fails the run.
#
# This is for local runs, where a mismatch against a reference generated on
# another machine is not worth failing over. Defaults to false: CI does not
# set it, so a mismatch there is a failure.
tolerate_hash_mismatch: bool = False

# Names of the tests whose mismatches were tolerated this session, reported by
# the `pytest_sessionfinish` hook so a green run still says plainly what was
# let through.
tolerated_mismatches: list[str] = []

# Test name -> component -> path (str) of the joined output file, for every
# test whose hashes differ from its reference, tolerated or not. Read by the
# `pytest_sessionfinish` hook, which writes it to mismatches_<environ>.json for
# the base-commit comparison in CI (tests/compare_outputs.py). A test with no
# reference entry has nothing to compare against and is not listed.
mismatched_outputs: dict[str, dict[str, str]] = {}

# Equal statistics with a changed hash mean either a change below about 1e-10
# relative or values that moved without changing min, max, mean or rms (a
# shifted or reordered field); only the per-point base comparison tells which.
STATISTICS_EQUAL = "statistics equal to 10 significant digits (see the base-commit report)"

UPDATE_HINT = (
    "To accept these results, run: python tests/update_references.py --run <RUN_ID>  "
    "(RUN_ID of the CI run whose artifacts hold the new values; see tests/README.md)"
)


def _sha256s(test_entry: dict) -> dict[str, dict[str, str]]:
    """``{component: {variable: SHA-256}}`` of a test entry, without the statistics."""
    return {
        component: {name: entry["sha256"] for name, entry in variables.items()}
        for component, variables in test_entry.items()
    }


def _describe_component_diff(name: str, expected: dict, computed: dict) -> list[str]:
    """Lines naming the variables whose hashes differ in one component, and how far they moved."""

    def show(entry: dict, stat: str) -> str:
        return f"{entry[stat]:.10g}" if stat in entry else "none"

    lines = [f"  {name}:"]
    for var, new in computed.items():
        if var not in expected or expected[var]["sha256"] == new["sha256"]:
            continue
        old = expected[var]
        # A statistic missing on one side (e.g. the variable became all-NaN) has moved too.
        stats = [stat for stat in dict.fromkeys([*old, *new]) if stat != "sha256"]
        moved = [
            f"{stat} {show(old, stat)} -> {show(new, stat)}"
            for stat in stats
            if stat not in old or stat not in new or old[stat] != new[stat]
        ]
        lines.append(
            f"    {var}: " + (", ".join(moved) or STATISTICS_EQUAL)
        )
    only_computed = [v for v in computed if v not in expected]
    only_reference = [v for v in expected if v not in computed]
    if only_computed:
        lines.append(f"    only in computed: {', '.join(only_computed)}")
    if only_reference:
        lines.append(f"    only in reference: {', '.join(only_reference)}")
    return lines


def assert_output_matches_reference(
    reference_results: dict, test_name: str, summaries: "dict[str, OutputSummary]"
):
    """Compare ``summaries`` against ``reference_results[test_name]``.

    ``summaries`` maps each output component to its :class:`OutputSummary`, e.g.
    ``{"physics": pv}`` or ``{"physics": pv, "bgc": bv}``. Only the SHA-256 of
    each variable is compared; the statistics next to it explain a mismatch.

    The computed values are recorded into :data:`computed_results` before the
    comparison is made, so the session-end dump captures every test that
    ran (not just the ones that failed).

    When :data:`tolerate_hash_mismatch` is set, a mismatch calls
    ``pytest.xfail`` rather than raising, which ends the test immediately and
    reports it as XFAIL. As with the ``AssertionError`` below, nothing after
    this call in the test body runs, so a test is only ever recorded once.
    """
    computed = {component: s.variables for component, s in summaries.items()}
    computed_results[test_name] = computed

    if test_name not in reference_results:
        message = (
            f"No reference results for {test_name!r} in results_{environ}.json.\n"
            f"{UPDATE_HINT}"
        )
    else:
        expected = reference_results[test_name]
        if any(predates_statistics(variables) for variables in expected.values()):
            # The hashes cannot be read per variable, so every output is suspect.
            differing = list(summaries)
            lines = [
                "  reference predates per-variable statistics; "
                "regenerate it with the command below"
            ]
        else:
            old, new = _sha256s(expected), _sha256s(computed)
            differing = [c for c in dict.fromkeys([*old, *new]) if old.get(c, {}) != new.get(c, {})]
            if not differing:
                return
            lines = [
                line
                for c in differing
                for line in _describe_component_diff(c, expected.get(c, {}), computed.get(c, {}))
            ]
        # A component present only in the reference has no head file to compare.
        paths = {c: str(summaries[c].path) for c in differing if c in summaries}
        if paths:
            mismatched_outputs[test_name] = paths
        message = (
            f"Hash mismatch for {test_name!r} (reference: results_{environ}.json).\n"
            + "\n".join(lines)
            + f"\n{UPDATE_HINT}"
        )

    if tolerate_hash_mismatch:
        tolerated_mismatches.append(test_name)
        pytest.xfail(message)
    raise AssertionError(message)
