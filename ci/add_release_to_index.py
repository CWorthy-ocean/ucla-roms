#!/usr/bin/env python3
"""
Add a new release to the release-notes index, ``docs/releases.rst``.

Each release has its own notes file, ``docs/releases/<tag>.rst``, and
``docs/releases.rst`` pulls them into one page with ``.. include::`` lines,
newest first. ``towncrier build`` writes the new file but knows nothing of the
index, so ``finalize_release_notes.yaml`` runs this script first.
``release_body.py --check-only`` later requires the release to be the top
include.

The script stops with an error if ``docs/releases/<tag>.rst`` already exists or
the index already includes it: towncrier would overwrite that file without
asking.

Release tags have no leading ``v`` (``0.15.0``); one on the argument is
stripped.

Usage:
    python ci/add_release_to_index.py <tag>
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
_INCLUDE_RE = re.compile(r"^\.\. include:: releases/\S+\.rst[ \t]*$", re.MULTILINE)


def normalize_tag(tag: str) -> str:
    """Strip a leading ``v``/``V`` from *tag* if present (``v0.7.0`` -> ``0.7.0``)."""
    return tag[1:] if re.match(r"^[vV]\d", tag) else tag


def add_include(index_text: str, tag: str) -> str:
    """
    Return *index_text* with an include of ``releases/<tag>.rst`` as the newest.

    The line goes above the first existing include, or at the end of an index
    that has none. Exits with an error if the index already includes *tag*.
    """
    line = f".. include:: releases/{tag}.rst"
    if line in index_text.splitlines():
        sys.exit(f"The release-notes index already includes releases/{tag}.rst.")
    first = _INCLUDE_RE.search(index_text)
    if first is None:
        return f"{index_text.rstrip()}\n\n{line}\n"
    return f"{index_text[: first.start()]}{line}\n\n{index_text[first.start() :]}"


def add_release(repo: Path, tag: str) -> None:
    """
    Add *tag* to ``docs/releases.rst`` under *repo*.

    Exits with an error, changing nothing, if the release already has a notes
    file or an include.
    """
    notes = repo / "docs" / "releases" / f"{tag}.rst"
    if notes.exists():
        sys.exit(
            f"{notes.relative_to(repo)} already exists; "
            "towncrier would overwrite it. Is this tag already released?"
        )
    index = repo / "docs" / "releases.rst"
    index.write_text(add_include(index.read_text(), tag))
    print(f"Added releases/{tag}.rst to {index.relative_to(repo)}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("tag", help="Release tag, e.g. 0.15.0")
    args = parser.parse_args()
    add_release(REPO_ROOT, normalize_tag(args.tag))


if __name__ == "__main__":
    main()
