"""Helpers for asserting test output hashes match the reference JSON.

Every call to :func:`assert_matches_reference` records its computed values
into the module-level :data:`computed_results` dict, regardless of whether
the assertion passed or failed. ``conftest.py``'s ``pytest_sessionfinish``
hook then dumps the collected dict at the end of the run, so you can paste
it straight into ``results/results_<environ>.json`` to update references.
"""

# Session-wide collector, populated by `assert_matches_reference` and read
# by the `pytest_sessionfinish` hook in conftest.py. Only tests that actually
# ran will appear here, so running a subset of the suite produces a partial
# dict rather than entries with missing values.
computed_results: dict[str, dict] = {}


def assert_output_matches_reference(reference_results: dict, test_name: str, computed: dict):
    """Compare ``computed`` against ``reference_results[test_name]``.

    ``computed`` is a dict like ``{"physics": "<hash>"}`` or
    ``{"physics": "<hash>", "bgc": "<hash>"}``.

    The computed dict is recorded into :data:`computed_results` before the
    comparison is made, so the session-end dump captures every test that
    ran (not just the ones that failed).
    """
    computed_results[test_name] = computed

    expected = reference_results.get(test_name)
    if expected != computed:
        raise AssertionError(
            f"Hash mismatch for {test_name!r}.\n"
            f"  expected: {expected}\n"
            f"  computed: {computed}"
        )
