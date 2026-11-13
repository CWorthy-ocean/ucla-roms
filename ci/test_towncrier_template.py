"""
Tests for the towncrier setup that turns fragments into ``docs/releases/<tag>.rst``.

Run with ``pytest ci/test_towncrier_template.py``. Needs ``towncrier`` and
``docutils`` and is skipped without them; like the other ``ci/`` tests, CI does
not run it.

Each test builds a release from sample fragments with the repository's own
``towncrier.toml`` and template, in a scratch git repository, exactly as
``finalize_release_notes.yaml`` does (index first, then ``towncrier build``).
It then checks the result against the style of the releases already in
``docs/releases/``, parses it as RST as ``changelog.yaml`` does, and reads it
back with ``release_body.py``, as ``publish_release.yaml`` does.
"""

from __future__ import annotations

import importlib.util
import io
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path
from types import ModuleType

import pytest

pytest.importorskip("towncrier")
docutils_core = pytest.importorskip("docutils.core")

REPO = Path(__file__).resolve().parent.parent
# The fragments folder, relative to the repository root (docs/changes).
CHANGES = tomllib.loads((REPO / "towncrier.toml").read_text())["tool"]["towncrier"][
    "directory"
]
URL = "https://github.com/CWorthy-ocean/ucla-roms/pull"
INDEX = "Release notes\n#############\n\n.. include:: releases/0.14.0.rst\n"
OLD_RELEASE = (
    ".. _0.14.0:\n\n0.14.0\n------\n\n"
    f"Bug Fixes\n~~~~~~~~~\n\n- An old fix. (`#1 <{URL}/1>`_)\n"
)
HEAD = ".. _0.15.0:\n\n0.15.0\n------\n\n"


def _load(name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, REPO / "ci" / f"{name}.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


release_body = _load("release_body")
add_release_to_index = _load("add_release_to_index")


def _link(pr: int) -> str:
    return f"(`#{pr} <{URL}/{pr}>`_)"


def _run(cwd: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(args, cwd=cwd, capture_output=True, text=True, check=False)


@pytest.fixture
def scratch(tmp_path: Path, monkeypatch) -> Path:
    """A git repository holding the real towncrier config and template."""
    shutil.copy(REPO / "towncrier.toml", tmp_path)
    (tmp_path / "ci").mkdir()
    shutil.copy(REPO / "ci" / "towncrier_template.rst", tmp_path / "ci")
    (tmp_path / "docs" / "releases").mkdir(parents=True)
    (tmp_path / "docs" / "releases.rst").write_text(INDEX)
    (tmp_path / "docs" / "releases" / "0.14.0.rst").write_text(OLD_RELEASE)
    (tmp_path / CHANGES).mkdir()
    shutil.copy(REPO / CHANGES / "README.md", tmp_path / CHANGES)
    _run(tmp_path, "git", "init", "-q")
    _run(tmp_path, "git", "add", "-A")
    monkeypatch.setattr(release_body, "RELEASES_DIR", tmp_path / "docs" / "releases")
    monkeypatch.setattr(release_body, "RELEASES_INDEX", tmp_path / "docs/releases.rst")
    return tmp_path


def _add_fragments(repo: Path, fragments: dict[str, str]) -> None:
    for name, text in fragments.items():
        (repo / CHANGES / name).write_text(text)
    _run(repo, "git", "add", "-A")


def _build(repo: Path, version: str) -> subprocess.CompletedProcess[str]:
    return _run(
        repo, sys.executable, "-m", "towncrier", "build", "--yes", "--version", version
    )


def _release(repo: Path, version: str = "0.15.0") -> str:
    """Finalize *version* as the workflow does and return its notes file."""
    add_release_to_index.add_release(repo, version)
    result = _build(repo, version)
    assert result.returncode == 0, result.stderr
    return (repo / "docs" / "releases" / f"{version}.rst").read_text()


def _assert_clean_rst(text: str) -> None:
    """Fail on anything docutils reports at warning level or above."""
    warnings = io.StringIO()
    docutils_core.publish_doctree(
        text, settings_overrides={"report_level": 2, "warning_stream": warnings}
    )
    assert warnings.getvalue() == ""


def test_build_matches_release_notes_style(scratch):
    _add_fragments(
        scratch,
        {
            "12.breaking.rst": "Lead ``text``.\n- First detail.\n- Second detail.\n",
            "16.output.rst": "Different numbers.\n",
            "14.deprecation.rst": "Deprecated: ``a``, removed in 1.0.0; use ``b``.\n",
            "12.feature.rst": "A feature.\n",
            # Not yet numbered: renders without a link until the release.
            "+later.feature.rst": "An unnumbered feature.\n",
            "13.bugfix.rst": "A later fix.\n",
            "11.bugfix.rst": "An earlier fix.\n- Its detail.\n",
            "15.misc.1.rst": "One note.\n",
            "15.misc.2.rst": "Another note.\n",
        },
    )
    text = _release(scratch)

    body = (
        "Breaking Changes\n~~~~~~~~~~~~~~~~\n\n"
        f"- Lead ``text``. {_link(12)}\n\n"
        "  - First detail.\n"
        "  - Second detail.\n\n"
        "Changes to Output\n~~~~~~~~~~~~~~~~~\n\n"
        f"- Different numbers. {_link(16)}\n\n"
        "Deprecations\n~~~~~~~~~~~~\n\n"
        f"- Deprecated: ``a``, removed in 1.0.0; use ``b``. {_link(14)}\n\n"
        "New Features\n~~~~~~~~~~~~\n\n"
        f"- A feature. {_link(12)}\n"
        "- An unnumbered feature.\n\n"
        "Bug Fixes\n~~~~~~~~~\n\n"
        f"- An earlier fix. {_link(11)}\n\n"
        "  - Its detail.\n\n"
        f"- A later fix. {_link(13)}\n\n"
        "Miscellaneous\n~~~~~~~~~~~~~\n\n"
    )
    assert text.startswith(HEAD + body)
    # Two notes from one PR and type come out in directory-listing order.
    misc = text[len(HEAD + body) :]
    assert sorted(misc.splitlines()) == [
        f"- Another note. {_link(15)}",
        f"- One note. {_link(15)}",
    ]
    _assert_clean_rst(text)
    assert [p.name for p in (scratch / CHANGES).iterdir()] == ["README.md"]

    # The new file is the top include, and the old release is untouched.
    assert (scratch / "docs" / "releases.rst").read_text() == INDEX.replace(
        ".. include:: releases/0.14.0.rst",
        ".. include:: releases/0.15.0.rst\n\n.. include:: releases/0.14.0.rst",
    )
    assert (scratch / "docs" / "releases" / "0.14.0.rst").read_text() == OLD_RELEASE

    # publish_release.yaml: the newest-release check, then the Markdown body.
    release_body.check_newest("0.15.0")
    markdown = release_body.section_body("0.15.0")
    assert markdown.startswith(
        "### Breaking Changes\n\n"
        f"- Lead `text`. ([#12]({URL}/12))\n\n"
        "  - First detail.\n"
        "  - Second detail.\n\n"
        "### Changes to Output\n\n"
    )
    assert "0.15.0" not in markdown and "~~" not in markdown


def test_fragment_written_as_a_bullet(scratch):
    _add_fragments(scratch, {"16.feature.rst": "- A feature.\n  - A detail.\n"})
    text = _release(scratch)
    assert text == (
        f"{HEAD}New Features\n~~~~~~~~~~~~\n\n"
        f"- A feature. {_link(16)}\n\n  - A detail.\n"
    )
    _assert_clean_rst(text)


def test_link_stays_on_the_first_paragraph(scratch):
    _add_fragments(
        scratch,
        {
            "17.feature.rst": "A feature.\n\n- A detail.\n",
            "18.feature.rst": "Wrapped\nlead.\n\nA second paragraph.\n",
            "19.bugfix.rst": "A fix, as in::\n\n    x = 1\n",
            "20.bugfix.rst": "A last fix.\n",
        },
    )
    text = _release(scratch)
    # One blank line before the sub-bullets, however many the fragment had.
    assert f"- A feature. {_link(17)}\n\n  - A detail.\n\n" in text
    assert f"- Wrapped\n  lead. {_link(18)}\n\n  A second paragraph.\n\n" in text
    # The link does not land in the literal block or on its "::" marker.
    assert (
        f"- A fix, as in::\n\n      x = 1\n\n  {_link(19)}\n\n- A last fix. {_link(20)}\n"
    ) in text
    assert not any(line != line.rstrip() for line in text.splitlines())
    _assert_clean_rst(text)


def test_build_without_fragments_is_not_empty(scratch):
    text = _release(scratch)
    assert text == f"{HEAD}No significant changes.\n"
    _assert_clean_rst(text)
    # release_body.py refuses an empty body, which would stop the release.
    assert release_body.section_body("0.15.0") == "No significant changes.\n"


@pytest.mark.parametrize("name", ["PR.misc.rst", "12.features.rst"])
def test_fragment_names_are_checked(scratch, name):
    _add_fragments(scratch, {name: "A note.\n"})
    result = _build(scratch, "0.15.0")
    assert result.returncode != 0
    assert not (scratch / "docs" / "releases" / "0.15.0.rst").exists()


def test_released_tag_is_refused_before_towncrier_overwrites_it(scratch):
    _add_fragments(scratch, {"12.feature.rst": "A feature.\n"})
    with pytest.raises(SystemExit, match=r"docs/releases/0\.14\.0\.rst already exists"):
        add_release_to_index.add_release(scratch, "0.14.0")
    assert (scratch / "docs" / "releases.rst").read_text() == INDEX
    assert (scratch / "docs" / "releases" / "0.14.0.rst").read_text() == OLD_RELEASE


@pytest.mark.parametrize(
    "note",
    ["The ub_ coefficients.\n", "Lead.\n- A detail.\nText that is not indented.\n"],
)
def test_invalid_rst_in_a_fragment_is_detectable(scratch, note):
    """What changelog.yaml's last step relies on: such a note builds, but warns."""
    _add_fragments(scratch, {"12.bugfix.rst": note})
    text = _release(scratch)
    with pytest.raises(AssertionError):
        _assert_clean_rst(text)
