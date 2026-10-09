"""Merge computed test results into the reference hash files in ``tests/results/``.

The reference hashes can only come from CI: the suite's inputs are generated on
the runner and the compilers differ from a laptop's, so a local run reproduces
neither. Each CI job uploads ``computed_<environ>.json`` as an artifact
(``results-github_gnu`` / ``results-github_ifx``). This script downloads those
artifacts with ``gh``, merges them into ``results/results_<environ>.json`` and
prints what changed, so the diff can be reviewed and committed with the PR.

    python tests/update_references.py --run <run-id>
    python tests/update_references.py --from-file computed_github_gnu.json --environ github_gnu

Merge rule: keys present in the computed dump replace the reference entry, all
other reference keys are kept (``--prune`` deletes them). Standard library only.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = REPO_ROOT / "tests" / "results"
ARTIFACTS = ("results-github_gnu", "results-github_ifx")
COMPUTED_PREFIX = "computed_"


def download_artifacts(run_id: str, repo: str, dest: Path) -> list[str]:
    """Download the results artifacts of a CI run into ``dest``; return the names found.

    One ``gh`` call per artifact, each into its own directory, so a missing
    compiler artifact is reported by name and never blocks the other, and a
    partially populated directory is never retried into.
    """
    if shutil.which("gh") is None:
        sys.exit("error: the GitHub CLI 'gh' was not found on PATH")

    found = []
    for name in ARTIFACTS:
        sub = dest / name
        sub.mkdir()
        cmd = ["gh", "run", "download", run_id, "-R", repo, "-n", name, "-D", str(sub)]
        if subprocess.run(cmd).returncode == 0:
            found.append(name)
        else:
            print(f"warning: artifact {name} not downloaded from run {run_id}", file=sys.stderr)
    if not found:
        sys.exit(f"error: no results artifact could be downloaded from run {run_id} in {repo}")
    return found


def find_computed(root: Path) -> dict[str, Path]:
    """Map environ -> computed file for every ``computed_<environ>.json`` under ``root``."""
    found: dict[str, Path] = {}
    for path in sorted(root.rglob(f"{COMPUTED_PREFIX}*.json")):
        environ = path.stem[len(COMPUTED_PREFIX):]
        if environ in found:
            sys.exit(f"error: two computed files for {environ}: {found[environ]} and {path}")
        found[environ] = path
    return found


def describe_changes(old: dict | str, new: dict | str) -> list[str]:
    """One line per component of a test whose reference differs from the computed value."""
    # A reference or computed test entry that is not a {component: ...} dict has no
    # per-component structure to describe.
    if not isinstance(old, dict) or not isinstance(new, dict):
        return ["format changed"]
    lines = []
    for comp in [*old, *(c for c in new if c not in old)]:
        if comp not in new:
            lines.append(f"{comp}: component removed")
        elif comp not in old:
            lines.append(f"{comp}: component added")
        elif old[comp] == new[comp]:
            continue
        elif isinstance(old[comp], str) or isinstance(new[comp], str):
            # Legacy references hold one digest per component, not per variable.
            lines.append(f"{comp}: format changed")
        else:
            o, n = old[comp], new[comp]
            parts = []
            changed = sorted(v for v in o.keys() & n.keys() if o[v] != n[v])
            only_new = sorted(n.keys() - o.keys())
            only_old = sorted(o.keys() - n.keys())
            if changed:
                parts.append(f"changed: {', '.join(changed)}")
            if only_new:
                parts.append(f"only in new: {', '.join(only_new)}")
            if only_old:
                parts.append(f"only in old: {', '.join(only_old)}")
            lines.append(f"{comp}: " + "; ".join(parts))
    return lines


def merge(environ: str, computed: dict, prune: bool, dry_run: bool) -> None:
    """Merge ``computed`` into the reference file for ``environ`` and print the report."""
    ref_path = RESULTS_DIR / f"results_{environ}.json"
    old = json.loads(ref_path.read_text()) if ref_path.exists() else {}

    # Existing keys keep their position; new keys are appended. Unsorted, so the
    # diff of the JSON file shows only what changed.
    merged = dict(old)
    merged.update(computed)
    removed = [k for k in old if k not in computed] if prune else []
    for key in removed:
        del merged[key]

    added = [k for k in computed if k not in old]
    changed = [k for k in computed if k in old and old[k] != computed[k]]
    unchanged = [k for k in computed if k in old and old[k] == computed[k]]

    print(f"== {environ} ({ref_path.relative_to(REPO_ROOT)}) ==")
    for key in added:
        print(f"  added   {key}")
    for key in changed:
        print(f"  changed {key}")
        for line in describe_changes(old[key], computed[key]):
            print(f"    {line}")
    for key in removed:
        print(f"  removed {key}")
    print(
        f"{environ}: {len(added)} added, {len(changed)} changed, "
        f"{len(removed)} removed, {len(unchanged)} unchanged"
    )

    if dry_run:
        print(f"{environ}: dry run, nothing written")
    elif merged != old or not ref_path.exists():
        ref_path.write_text(json.dumps(merged, indent=2) + "\n")
        print(f"{environ}: wrote {ref_path.relative_to(REPO_ROOT)}")
    else:
        print(f"{environ}: reference file already up to date")


def load_computed(path: Path) -> dict:
    """Read a computed dump, checking it has the ``{test: {component: ...}}`` shape."""
    data = json.loads(path.read_text())
    if not isinstance(data, dict) or not all(isinstance(v, dict) for v in data.values()):
        sys.exit(f"error: {path} is not a {{test: {{component: ...}}}} mapping")
    return data


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge computed test results from CI into tests/results/results_<environ>.json."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--run", metavar="RUN_ID", help="GitHub Actions run id to download the results artifacts from")
    source.add_argument("--from-file", metavar="PATH", type=Path, help="merge this local computed_<environ>.json instead")
    parser.add_argument("--repo", default="CWorthy-ocean/ucla-roms", help="repository of --run (default: %(default)s)")
    parser.add_argument("--environ", help="environment name for --from-file, e.g. github_gnu")
    parser.add_argument(
        "--prune",
        action="store_true",
        help="also delete reference keys absent from the computed dump; only meaningful "
        "after a full-suite run, since a -k subset would otherwise delete valid references",
    )
    parser.add_argument("--dry-run", action="store_true", help="print the report, write nothing")
    args = parser.parse_args()

    if args.from_file is not None:
        if not args.environ:
            parser.error("--from-file requires --environ")
        computed = {args.environ: load_computed(args.from_file)}
    else:
        if args.environ:
            parser.error("--environ is only used with --from-file")
        with tempfile.TemporaryDirectory() as tmp:
            found = download_artifacts(args.run, args.repo, Path(tmp))
            print(f"downloaded: {', '.join(found)}")
            files = find_computed(Path(tmp))
            if not files:
                sys.exit(f"error: no {COMPUTED_PREFIX}<environ>.json found in the artifacts of run {args.run}")
            computed = {environ: load_computed(path) for environ, path in files.items()}

    for environ, data in computed.items():
        merge(environ, data, args.prune, args.dry_run)


if __name__ == "__main__":
    main()
