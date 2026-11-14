"""
Tests for ``add_release_to_index.py`` (a new release's include in the index).

Run with ``pytest ci/test_add_release_to_index.py``; like the other ``ci/``
tests, CI does not run it.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent

_SPEC = importlib.util.spec_from_file_location(
    "add_release_to_index", REPO / "ci" / "add_release_to_index.py"
)
assert _SPEC is not None and _SPEC.loader is not None
ari = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(ari)

_INDEX = (
    "Release notes\n#############\n\n"
    ".. include:: releases/0.14.0.rst\n\n"
    ".. include:: releases/0.13.1.rst"
)


def test_new_release_goes_above_the_others():
    assert ari.add_include(_INDEX, "0.15.0") == (
        "Release notes\n#############\n\n"
        ".. include:: releases/0.15.0.rst\n\n"
        ".. include:: releases/0.14.0.rst\n\n"
        ".. include:: releases/0.13.1.rst"
    )


def test_index_without_includes():
    assert ari.add_include("Release notes\n#############\n", "0.1.0") == (
        "Release notes\n#############\n\n.. include:: releases/0.1.0.rst\n"
    )


@pytest.mark.parametrize("tag", ["0.14.0", "0.13.1"])
def test_included_release_is_refused(tag):
    with pytest.raises(SystemExit, match=f"already includes releases/{tag}.rst"):
        ari.add_include(_INDEX, tag)


def test_a_longer_tag_is_not_mistaken_for_an_included_one():
    assert ".. include:: releases/0.14.0.1.rst\n\n" in ari.add_include(
        _INDEX, "0.14.0.1"
    )


def test_normalize_tag():
    assert ari.normalize_tag("v0.15.0") == "0.15.0"
    assert ari.normalize_tag("0.15.0") == "0.15.0"


def test_add_release_writes_the_index(tmp_path):
    (tmp_path / "docs" / "releases").mkdir(parents=True)
    index = tmp_path / "docs" / "releases.rst"
    index.write_text(_INDEX)
    ari.add_release(tmp_path, "0.15.0")
    assert index.read_text().splitlines()[3] == ".. include:: releases/0.15.0.rst"


def test_existing_notes_file_is_refused(tmp_path):
    (tmp_path / "docs" / "releases").mkdir(parents=True)
    index = tmp_path / "docs" / "releases.rst"
    index.write_text("Release notes\n#############\n")
    # A notes file the index does not include: towncrier would still overwrite it.
    (tmp_path / "docs" / "releases" / "0.15.0.rst").write_text("0.15.0\n------\n")
    with pytest.raises(SystemExit, match=r"0\.15\.0\.rst already exists"):
        ari.add_release(tmp_path, "0.15.0")
    assert index.read_text() == "Release notes\n#############\n"


def test_real_index_takes_a_new_release():
    """The repository's own index has the shape ``add_include`` expects."""
    text = (REPO / "docs" / "releases.rst").read_text()
    new = ari.add_include(text, "999.0.0")
    assert new.replace(".. include:: releases/999.0.0.rst\n\n", "", 1) == text
    first = next(ln for ln in new.splitlines() if ln.startswith(".. include::"))
    assert first == ".. include:: releases/999.0.0.rst"
