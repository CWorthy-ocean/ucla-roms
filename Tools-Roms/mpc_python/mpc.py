#!/usr/bin/env python3
"""
Python re-implementation of the MPC (multi-functional pre-compiler)
Originally by Alexander Shchepetkin
Rewritten in Python with modular passes and switchable features.
"""

import sys
import argparse
from pathlib import Path

# Import transformation passes
from passes.cleanup import pass_cleanup
from passes.trap_unmatched_quotes import pass_trap_unmatched_quotes
from passes.f77_to_f95 import pass_f77_to_f95
from passes.double_constants import pass_double_constants
from passes.zig_zag import pass_zig_zag
from passes.trap_barriers import pass_trap_barriers
from passes.fold_lines import pass_fold_lines
from passes.emit import pass_emit

class MPCConfig:
    """Holds all feature switches."""
    def __init__(self,
                 trap_unmatched_quotes=True,
                 f77_to_f95=True,
                 double_constants=True,
                 zig_zag=True,
                 trap_barriers=False,
                 fold_lines=True):
        self.trap_unmatched_quotes = trap_unmatched_quotes
        self.f77_to_f95 = f77_to_f95
        self.double_constants = double_constants
        self.zig_zag = zig_zag
        self.trap_barriers = trap_barriers
        self.fold_lines = fold_lines


def build_pipeline(cfg: MPCConfig):
    """Return an ordered list of transformation passes."""
    pipeline = []

    pipeline.append(pass_cleanup)

    if cfg.trap_unmatched_quotes:
        pipeline.append(pass_trap_unmatched_quotes)

    if cfg.f77_to_f95:
        pipeline.append(pass_f77_to_f95)

    if cfg.double_constants:
        pipeline.append(pass_double_constants)

    if cfg.zig_zag:
        pipeline.append(pass_zig_zag)

    if cfg.trap_barriers:
        pipeline.append(pass_trap_barriers)

    if cfg.fold_lines:
        pipeline.append(pass_fold_lines)

    pipeline.append(pass_emit)

    return pipeline


def run_mpc(lines, cfg: MPCConfig):
    """Run all enabled passes on the list of source lines."""
    pipeline = build_pipeline(cfg)

    for transform in pipeline:
        lines = transform(lines, cfg)

    return lines


def main():
    parser = argparse.ArgumentParser(description="Python MPC pre-compiler")
    parser.add_argument("infile", nargs="?", default="-")
    parser.add_argument("outfile", nargs="?", default="-")

    # Feature switches
    parser.add_argument("--trap-unmatched-quotes", action="store_false")
    parser.add_argument("--f77-to-f95", action="store_false")
    parser.add_argument("--double-constants", action="store_false")
    parser.add_argument("--zig-zag", action="store_false")
    parser.add_argument("--trap-barriers", action="store_true")
    parser.add_argument("--no-fold", action="store_true")

    args = parser.parse_args()

    cfg = MPCConfig(
        trap_unmatched_quotes=args.trap_unmatched_quotes,
        f77_to_f95=args.f77_to_f95,
        double_constants=args.double_constants,
        zig_zag=args.zig_zag,
        trap_barriers=args.trap_barriers,
        fold_lines=not args.no_fold
    )

    # Input
    if args.infile == "-":
        lines = sys.stdin.read().splitlines()
    else:
        lines = Path(args.infile).read_text().splitlines()

    result = run_mpc(lines, cfg)

    # Output
    if args.outfile == "-":
        for L in result:
            print(L)
    else:
        Path(args.outfile).write_text("\n".join(result))


if __name__ == "__main__":
    main()
