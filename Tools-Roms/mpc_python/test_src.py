
# conftest.py
import os
import glob
from mpc import MPCConfig, run_mpc
import pytest
import tempfile
import shutil
from pathlib import Path
import subprocess


# Commit from which to fetch legacy MPC code for comparison
COMMIT = "75208625d3e5afd96fa9db8e809146da4c53d7de"
roms_root = os.environ.get("ROMS_ROOT","")

@pytest.fixture(scope="session")
def legacy_mpc_binary(tmp_path_factory):
    """
    Build the historical MPC exactly once per test session.
    Must use tmp_path_factory because session-scoped fixtures
    cannot depend on function-scoped fixtures like tmp_path.
    """

    build_dir = tmp_path_factory.mktemp("legacy_mpc")

    def extract(commit_path, local_name):
        result = subprocess.run(
            ["git", "show", f"{COMMIT}:{commit_path}"],
            text=True,
            capture_output=True,
            check=True,
        )
        out = build_dir / local_name
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(result.stdout)
        return out

    # Historical source layout
    extract("Tools-Roms/mpc.F",       "mpc.F")
    extract("Tools-Roms/Makefile",    "Makefile")
    extract("Tools-Roms/Make.depend", "Make.depend")
    extract("src/Makedefs.inc",       "../src/Makedefs.inc")

    # Build once
    subprocess.run(
        ["make", "mpc"],
        cwd=build_dir,
        check=True,
        text=True,
        capture_output=True
    )

    return str(build_dir / "mpc")

# Collect all .F test inputs
# test_inputs = sorted(glob.glob("../src//*.F"))
test_inputs = sorted(
    glob.glob(roms_root + "/src/*.F") +
    glob.glob(roms_root + "/src/*.h")
)

@pytest.mark.parametrize("fname", test_inputs)
def test_src(fname, legacy_mpc_binary,tmp_path):

    fname=Path(fname)
    fname_base = Path(fname).name

    binary = legacy_mpc_binary
    old=subprocess.run(f"{binary} {fname_base}",
                       text=True, shell=True, check=True, capture_output=True,
                       cwd=roms_root+"/src")

    expected=old.stdout.splitlines(keepends=True)


    # Read input .F file
    file_in = Path(fname).read_text().splitlines()

    # Build MPC config
    cfg = MPCConfig(
        trap_unmatched_quotes=True,
        f77_to_f95=True,
        double_constants=True,
        zig_zag=True,
        trap_barriers=False,
        fold_lines=True,
    )

    # Run python MPC
    file_out = run_mpc(file_in, cfg)
    actual = [line + "\n" for line in file_out]

    # Determine expected output file name
    # expected = Path(expected_path).read_text().splitlines(keepends=True)

    # Assertion with nice error message
    assert actual == expected, (
        f"\nFAILED on {fname}\n"
        f"\nEXPECTED:\n{''.join(expected)}\n"
        f"\nGOT:\n{''.join(actual)}"
    )
