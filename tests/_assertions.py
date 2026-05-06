"""Helpers for asserting test output hashes match the reference JSON."""


def assert_output_matches_reference(reference_results: dict, test_name: str, computed: dict):
    """Compare ``computed`` against ``reference_results[test_name]``.

    ``computed`` is a dict like ``{"physics": "<hash>"}`` or
    ``{"physics": "<hash>", "bgc": "<hash>"}``.

    On mismatch the assertion message shows the full computed block as a
    JSON-paste-ready snippet so the reference file can be updated quickly
    when an intentional change lands.
    """
    expected = reference_results.get(test_name)
    if expected != computed:
        import json
        snippet = json.dumps({test_name: computed}, indent=2)
        raise AssertionError(
            f"Hash mismatch for {test_name!r}.\n"
            f"  expected: {expected}\n"
            f"  computed: {computed}\n"
            f"To update the reference file, paste:\n{snippet}"
        )
