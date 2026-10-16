"""Print what moved between two ROMS output files, variable by variable.

When a CI run fails on hash mismatches, a SHA-256 mismatch says only that a
variable changed. This script opens the head output and the same test's output
from the PR's base commit and reports, per variable, how many points differ,
the largest absolute difference, that difference relative to the largest base
value, and where it sits, so a reviewer can tell round-off from a real change
without opening NetCDF. It is a report, not a gate: it exits 0 whenever it
produced one.

    python tests/compare_outputs.py BASE.nc HEAD.nc
    python tests/compare_outputs.py --mismatches tests/results/mismatches_github_gnu.json \\
        --head-root /tmp/pytest_head --base-root /tmp/pytest_base [--markdown OUT]

The mismatches file (written by the test session when hashes disagree) maps
test -> component -> absolute path of the head output file. The head run and
the base run use different ``--basetemp`` directories and pytest names a test's
tmp directory after the test, so the base file is the head path with
``--head-root`` replaced by ``--base-root``. That relies on test names (and
fixture ``mktemp`` names) being unique in their first 30 characters: pytest
numbers repeated names in run order, which could pair different tests. ``--markdown`` appends a Markdown
version of the report to OUT, which CI points at ``$GITHUB_STEP_SUMMARY``.

Needs numpy and xarray only.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr

NO_BASE = "no base output (test added or renamed, or the base run failed)"
NOT_AVAILABLE = "n/a"

# Sort groups for a changed variable: structural changes come first because
# they say the two runs wrote different things, which a magnitude cannot.
# A variable that is not numeric has no magnitude, so it sorts with them.
GROUP_ONLY_ONE_FILE = 0
GROUP_SHAPE_OR_DTYPE = 1
GROUP_NONFINITE = 2
GROUP_NON_NUMERIC = 3
GROUP_NUMERIC = 4


@dataclass
class VarDiff:
    """One variable that differs between base and head."""

    name: str
    group: int
    points: str  # the "points differing" cell: a count, or the structural change
    max_abs: float | None = None  # None: no position where both sides are finite
    rel: float | None = None  # max_abs / max |base|; inf when the base is all zero
    location: str = ""


@dataclass
class Comparison:
    """All variables of a pair of files: the ones that differ, sorted, and how many are identical."""

    total: int
    changed: list[VarDiff]

    @property
    def identical(self) -> int:
        return self.total - len(self.changed)


class MissingOutput(Exception):
    """A head or base file the comparison needs is not there; the message says why."""


def format_shape(shape: tuple[int, ...]) -> str:
    return "(" + ", ".join(str(n) for n in shape) + ")"


def format_location(dims: tuple, index: tuple[int, ...]) -> str:
    """``time=1, s_rho=9, ...`` for a position; ``scalar`` for a variable with no dimensions."""
    return ", ".join(f"{d}={i}" for d, i in zip(dims, index)) or "scalar"


def compare_numeric(name: str, dims: tuple, base: np.ndarray, head: np.ndarray) -> VarDiff | None:
    """Compare two numeric arrays of equal shape in float64; None when they are identical."""
    base = np.asarray(base, dtype=np.float64)
    head = np.asarray(head, dtype=np.float64)

    # NaN != NaN, so a plain != would flag every never-written fill value (NaN
    # once mask_and_scale has run); NaN on both sides counts as equal.
    base_nan = np.isnan(base)
    head_nan = np.isnan(head)
    differs = (base != head) & ~(base_nan & head_nan)
    n_differ = int(differs.sum())
    if n_differ == 0:
        return None

    # The magnitude is judged only where both values are finite. A differing
    # position where either side is NaN or inf (a blow-up, or a fill value that
    # appeared or vanished) has no size; it is counted and sorts first.
    both_finite = np.isfinite(base) & np.isfinite(head)
    nonfinite = int((differs & ~both_finite).sum())
    points = f"{n_differ}/{base.size}"
    if nonfinite:
        points += f" ({nonfinite} NaN or inf)"

    with np.errstate(invalid="ignore"):
        delta = np.where(differs & both_finite, np.abs(head - base), -1.0)
    flat = int(delta.argmax())
    max_abs = float(delta.flat[flat])
    if max_abs < 0:
        return VarDiff(name, GROUP_NONFINITE, points)

    location = format_location(dims, np.unravel_index(flat, delta.shape))
    finite_base = np.isfinite(base)
    max_base = float(np.abs(base[finite_base]).max()) if finite_base.any() else 0.0
    # A base that is zero everywhere (a CDR tracer before its release) has no
    # scale to be relative to: infinite here sorts it first, and it prints n/a.
    rel = max_abs / max_base if max_base > 0 else float("inf")
    group = GROUP_NONFINITE if nonfinite else GROUP_NUMERIC
    return VarDiff(name, group, points, max_abs, rel, location)


def compare_variable(name: str, base: xr.DataArray | None, head: xr.DataArray | None) -> VarDiff | None:
    """Describe how one variable differs between base and head; None when it does not."""
    if base is None:
        return VarDiff(name, GROUP_ONLY_ONE_FILE, "only in head")
    if head is None:
        return VarDiff(name, GROUP_ONLY_ONE_FILE, "only in base")
    if base.shape != head.shape:
        return VarDiff(name, GROUP_SHAPE_OR_DTYPE, f"shape {format_shape(base.shape)} -> {format_shape(head.shape)}")

    base_numeric = base.dtype.kind in "iuf"
    head_numeric = head.dtype.kind in "iuf"
    if base_numeric != head_numeric:
        return VarDiff(name, GROUP_SHAPE_OR_DTYPE, f"dtype {base.dtype} -> {head.dtype}")

    if base_numeric:
        diff = compare_numeric(name, base.dims, base.values, head.values)
        # The suite hashes raw bytes, so a dtype change moves the hash even when
        # every value compares equal in float64; say so rather than report "identical".
        if diff is None and base.dtype != head.dtype:
            return VarDiff(name, GROUP_SHAPE_OR_DTYPE, f"dtype {base.dtype} -> {head.dtype} (values equal)")
        return diff

    if np.array_equal(base.values, head.values):
        return None
    return VarDiff(name, GROUP_NON_NUMERIC, "differs")


def sort_key(diff: VarDiff) -> tuple:
    """Structural groups first; numeric changes by relative size, then absolute size, both descending."""
    if diff.group != GROUP_NUMERIC:
        return (diff.group, 0.0, 0.0, diff.name)
    rel = float("inf") if diff.rel is None else diff.rel
    return (diff.group, -rel, -(diff.max_abs or 0.0), diff.name)


def compare_files(base_path: Path, head_path: Path) -> Comparison:
    """Compare every data variable of two NetCDF files."""
    # Raw time numbers are enough for a difference. Decoding would also fail on the
    # never-written avg_begin_time/avg_end_time, whose fill values carry time units.
    opts = dict(decode_times=False, decode_timedelta=False)
    with xr.open_dataset(base_path, **opts) as base, xr.open_dataset(head_path, **opts) as head:
        names = sorted(base.data_vars.keys() | head.data_vars.keys())
        diffs = [compare_variable(n, base.data_vars.get(n), head.data_vars.get(n)) for n in names]
    changed = sorted((d for d in diffs if d is not None), key=sort_key)
    return Comparison(total=len(names), changed=changed)


def pair_paths(recorded: str, head_root: Path, base_root: Path) -> tuple[Path, Path]:
    """Map a recorded head output path to ``(base, head)`` file paths, or raise MissingOutput."""
    # Resolve before relative_to: pytest resolves --basetemp, so on macOS a
    # recorded path is /private/tmp/... while a root typed as /tmp/... is not.
    head = Path(recorded).resolve()
    try:
        rel = head.relative_to(head_root)
    except ValueError:
        raise MissingOutput(f"head output {head} is not under the head root {head_root}") from None
    if not head.is_file():
        raise MissingOutput(f"head output file is missing: {head}")
    base = base_root / rel
    if not base.is_file():
        raise MissingOutput(NO_BASE)
    return base, head


def compare_mismatches(
    mismatches: dict[str, dict[str, str]], head_root: Path, base_root: Path
) -> list[tuple[str, Comparison | str]]:
    """Compare every recorded test and component; each entry is a Comparison or a message saying why not."""
    head_root = head_root.resolve()
    base_root = base_root.resolve()
    results: list[tuple[str, Comparison | str]] = []
    for test, components in mismatches.items():
        for component, recorded in components.items():
            label = f"{test} / {component}"
            try:
                base, head = pair_paths(recorded, head_root, base_root)
                results.append((label, compare_files(base, head)))
            except MissingOutput as err:
                results.append((label, str(err)))
            except Exception as err:  # one unreadable file must not hide the other reports
                reason = str(err).splitlines()[0] if str(err) else ""
                results.append((label, f"could not compare: {type(err).__name__}: {reason}"))
    return results


def summary_line(comparison: Comparison) -> str:
    if not comparison.changed:
        return f"all {comparison.total} variables identical"
    return f"{len(comparison.changed)} of {comparison.total} variables differ"


def format_number(value: float | None) -> str:
    """``%.3g``; ``n/a`` when there is no finite value (or no scale, for infinity)."""
    if value is None or not np.isfinite(value):
        return NOT_AVAILABLE
    return f"{value:.3g}"


def table_cells(diff: VarDiff) -> list[str]:
    """The five table cells of a changed variable; structural changes leave the numeric ones empty."""
    if diff.group in (GROUP_ONLY_ONE_FILE, GROUP_SHAPE_OR_DTYPE, GROUP_NON_NUMERIC):
        return [diff.name, diff.points, "", "", ""]
    return [diff.name, diff.points, format_number(diff.max_abs), format_number(diff.rel), diff.location or NOT_AVAILABLE]


def render_text(label: str, result: Comparison | str) -> list[str]:
    """Text report of one test/component (or of a single pair, with an empty label)."""
    prefix = f"{label}: " if label else ""
    if isinstance(result, str):
        return [f"{prefix}{result}"]
    lines = [f"{prefix}{summary_line(result)}"]
    if result.changed:
        header = ["variable", "points differing", "max |Δ|", "rel", "at"]
        rows = [header] + [table_cells(d) for d in result.changed]
        # The last column is free text of varying length and is not padded.
        widths = [max(len(row[i]) for row in rows) for i in range(4)]
        for row in rows:
            padded = [cell.ljust(w) for cell, w in zip(row, widths)]
            lines.append(("  " + "  ".join([*padded, row[4]])).rstrip())
        lines.append(f"  identical: {result.identical} variables")
    return lines


def render_markdown(mismatches_file: Path, results: list[tuple[str, Comparison | str]]) -> list[str]:
    """Markdown report for the GitHub step summary."""
    lines = [
        "## Output differences from the base commit",
        "",
        f"Head outputs listed in `{mismatches_file.name}` compared with the same tests' outputs on the base commit.",
        "",
    ]
    for label, result in results:
        lines += [f"### {label}", ""]
        if isinstance(result, str):
            lines += [result, ""]
            continue
        lines += [summary_line(result), ""]
        if result.changed:
            lines += [
                "| variable | points differing | max abs Δ | max abs Δ / max abs base | at |",
                "|---|---|---|---|---|",
            ]
            for diff in result.changed:
                name, points, max_abs, rel, location = table_cells(diff)
                at = f"`{location}`" if location not in ("", NOT_AVAILABLE) else location
                # A bytes dtype prints as |S1, and a pipe in a cell would split the row.
                cells = [c.replace("|", "\\|") for c in (name, points, max_abs, rel, at)]
                lines.append("| " + " | ".join(cells) + " |")
            lines += ["", f"identical: {result.identical} variables", ""]
    return lines


def load_mismatches(path: Path) -> dict[str, dict[str, str]]:
    """Read a mismatches file, checking it has the ``{test: {component: path}}`` shape."""
    try:
        data = json.loads(path.read_text())
    except (OSError, ValueError) as err:
        sys.exit(f"error: cannot read mismatches file {path}: {err}")
    ok = isinstance(data, dict) and all(
        isinstance(v, dict) and all(isinstance(p, str) for p in v.values()) for v in data.values()
    )
    if not ok:
        sys.exit(f"error: {path} is not a {{test: {{component: path}}}} mapping")
    return data


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Per-variable report of how a test's head output differs from its base-commit output.",
        usage="%(prog)s BASE.nc HEAD.nc\n"
        "       %(prog)s --mismatches FILE --head-root DIR --base-root DIR [--markdown OUT]",
    )
    parser.add_argument("files", nargs="*", metavar="NC", help="BASE.nc HEAD.nc: compare this one pair")
    parser.add_argument("--mismatches", metavar="FILE", type=Path, help="mismatches_<environ>.json written by the head test run")
    parser.add_argument("--head-root", metavar="DIR", type=Path, help="--basetemp of the head run (e.g. /tmp/pytest_head)")
    parser.add_argument("--base-root", metavar="DIR", type=Path, help="--basetemp of the base run (e.g. /tmp/pytest_base)")
    parser.add_argument("--markdown", metavar="OUT", type=Path, help="also append a Markdown report to OUT (e.g. $GITHUB_STEP_SUMMARY)")
    args = parser.parse_args()

    if args.mismatches is None:
        if len(args.files) != 2:
            parser.error("give BASE.nc HEAD.nc, or --mismatches FILE with --head-root and --base-root")
        if args.head_root or args.base_root or args.markdown:
            parser.error("--head-root, --base-root and --markdown are only used with --mismatches")
        base, head = (Path(f) for f in args.files)
        for path in (base, head):
            if not path.is_file():
                sys.exit(f"error: {path} is not a file")
        print(f"base: {base}")
        print(f"head: {head}")
        print("\n".join(render_text("", compare_files(base, head))))
        return

    if args.files:
        parser.error("BASE.nc HEAD.nc and --mismatches are mutually exclusive")
    if args.head_root is None or args.base_root is None:
        parser.error("--mismatches requires --head-root and --base-root")

    results = compare_mismatches(load_mismatches(args.mismatches), args.head_root, args.base_root)
    if not results:
        print(f"no mismatches recorded in {args.mismatches}")
    for i, (label, result) in enumerate(results):
        if i:
            print()
        print("\n".join(render_text(label, result)))
    if args.markdown is not None:
        with open(args.markdown, "a", encoding="utf-8") as out:
            # The leading blank line keeps the heading apart from whatever earlier steps wrote.
            out.write("\n" + "\n".join(render_markdown(args.mismatches, results)) + "\n")


if __name__ == "__main__":
    main()
