"""Unit test of src/cdf_direct.F90, the direct (O_DIRECT) reader of classic netCDF files that pio_roms
uses for the start-up files when pio_direct_read is on.

Writes small files in every classic format (CDF-1, CDF-2, CDF-5) with odd shapes, every converted type,
special values, a large header and record variables; compiles tests/additional_files/test_cdf_direct.F90
with the build's compiler settings (src/Makedefs.inc); and checks that every variable read with
cdf_direct is bit for bit what nf90_get_var returns. Also: a netCDF-4 file and a truncated file are
refused, and a file with a single record variable (packed records) is read right."""
import os
import shutil
import subprocess
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from ._helpers import MAKE_CMD

ROMS_ROOT = Path(os.environ["ROMS_ROOT"])
SRC = ROMS_ROOT / "src"
TEST_F90 = ROMS_ROOT / "tests" / "additional_files" / "test_cdf_direct.F90"

MAKEFILE = """include {src}/Makedefs.inc
test_cdf_direct: cdf_direct.F90 test_cdf_direct.F90
\t$(CPP) -P $(CPPFLAGS) cdf_direct.F90 > cdf_direct.f90
\t$(CPP) -P $(CPPFLAGS) test_cdf_direct.F90 > test_cdf_direct.f90
\t$(CFT) $(FFLAGS) -o test_cdf_direct cdf_direct.f90 test_cdf_direct.f90 $(LCDF)
test_cdf_direct_O3: test_cdf_direct
\t$(CFT) $(FFLAGS) -O3 -o test_cdf_direct_O3 cdf_direct.f90 test_cdf_direct.f90 $(LCDF)
"""


def _special_doubles(n):
    vals = np.array([np.nan, -0.0, 0.0, 5e-324, -2.2e-308, 1.7e308, -1.7e308, np.inf, -np.inf, 1.0 / 3.0])
    return np.resize(vals, n)


def _write(path, fmt, big_attr=False):
    rng = np.random.default_rng(7)
    nx, ny, nz = 37, 23, 5          # no length a multiple of 4096 bytes
    with netCDF4.Dataset(path, "w", format=fmt) as ds:
        if big_attr:
            ds.history = "x" * 400_000  # header longer than the first 256-KiB read
        ds.createDimension("time", None)
        ds.createDimension("s_rho", nz)
        ds.createDimension("eta_rho", ny)
        ds.createDimension("xi_rho", nx)
        ds.createDimension("auxil", 6)
        t = ds.createVariable("ocean_time", "f8", ("time",))
        temp = ds.createVariable("temp", "f8", ("time", "s_rho", "eta_rho", "xi_rho"))
        salt = ds.createVariable("salt", "f4", ("time", "s_rho", "eta_rho", "xi_rho"))
        zeta = ds.createVariable("zeta", "f8", ("time", "eta_rho", "xi_rho"))
        step = ds.createVariable("time_step", "i4", ("time", "auxil"))
        h = ds.createVariable("h", "f8", ("eta_rho", "xi_rho"))
        mask = ds.createVariable("mask", "i4", ("eta_rho", "xi_rho"))
        sh = ds.createVariable("sh", "i2", ("eta_rho", "xi_rho"))
        by = ds.createVariable("by", "i1", ("eta_rho", "xi_rho"))
        sh3 = ds.createVariable("sh3", "i2", ("time", "s_rho", "eta_rho", "xi_rho"))
        for v in (temp, salt, zeta, h, sh3):
            v.set_auto_maskandscale(False)
        for r in range(2):
            t[r] = 86400.0 * (r + 1)
            temp[r] = rng.normal(size=(nz, ny, nx))
            temp[r, 0, 0, :10] = _special_doubles(10)
            s = rng.normal(size=(nz, ny, nx)).astype("f4")
            s[0, 0, :4] = np.array([np.nan, -0.0, 1e-45, 3.4e38], dtype="f4")
            salt[r] = s
            zeta[r] = rng.normal(size=(ny, nx))
            step[r] = np.arange(6) + 10 * r
            sh3[r] = rng.integers(-32768, 32767, size=(nz, ny, nx), dtype="i2")
        h[:] = _special_doubles(ny * nx).reshape(ny, nx)
        mask[:] = rng.integers(-(2**31), 2**31 - 1, size=(ny, nx), dtype="i4")
        sh[:] = rng.integers(-32768, 32767, size=(ny, nx), dtype="i2")
        by[:] = rng.integers(-128, 127, size=(ny, nx), dtype="i1")
        if fmt == "NETCDF3_64BIT_DATA":
            us = ds.createVariable("us", "u2", ("eta_rho", "xi_rho"))   # not converted: left to netCDF
            us[:] = rng.integers(0, 65535, size=(ny, nx), dtype="u2")


def _write_one_record_var(path):
    # a single record variable of an odd byte length: netCDF packs the records without padding
    with netCDF4.Dataset(path, "w", format="NETCDF3_64BIT_DATA") as ds:
        ds.createDimension("time", None)
        ds.createDimension("n", 7)
        v = ds.createVariable("v", "i2", ("time", "n"))
        for r in range(3):
            v[r] = np.arange(7, dtype="i2") * (r + 1) - 3


@pytest.fixture(scope="module")
def build_dir(tmp_path_factory):
    # built with the suite's make command (the CI ifx job adds COMPILER=intel to it), at the test build's
    # flags and also at -O3 as ROMS runs are built (gfortran -O3 once miscompiled a type-punned conversion)
    d = tmp_path_factory.mktemp("cdf_direct_build")
    shutil.copy(SRC / "cdf_direct.F90", d)
    shutil.copy(TEST_F90, d)
    (d / "Makefile").write_text(MAKEFILE.format(src=SRC))
    res = subprocess.run(f"{MAKE_CMD} test_cdf_direct test_cdf_direct_O3", shell=True, cwd=d,
                         capture_output=True, text=True)
    assert res.returncode == 0, res.stdout + res.stderr
    return d


@pytest.fixture(params=["test", "O3"])
def test_binary(request, build_dir):
    return build_dir / ("test_cdf_direct" if request.param == "test" else "test_cdf_direct_O3")


def _run(binary, path, mode):
    args = [str(binary), str(path)] + (["buffered"] if mode == "buffered" else [])
    res = subprocess.run(args, capture_output=True, text=True)
    return res.returncode, res.stdout + res.stderr


def _mode(binary, path):
    # O_DIRECT where the file system has it (not tmpfs); the buffered path checks the same code otherwise
    rc, out = _run(binary, path, "direct")
    if "open with O_DIRECT failed" in out:
        return "buffered"
    return "direct"


@pytest.mark.parametrize("fmt", ["NETCDF3_CLASSIC", "NETCDF3_64BIT_OFFSET", "NETCDF3_64BIT_DATA"])
@pytest.mark.parametrize("big_attr", [False, True])
def test_classic_formats(test_binary, tmp_path, fmt, big_attr):
    path = tmp_path / f"t_{fmt}.nc"
    _write(path, fmt, big_attr)
    rc, out = _run(test_binary, path, _mode(test_binary, path))
    assert rc == 0 and "RESULT" in out and " failures 0 " in out, out


def test_packed_single_record_variable(test_binary, tmp_path):
    path = tmp_path / "one_rec.nc"
    _write_one_record_var(path)
    rc, out = _run(test_binary, path, _mode(test_binary, path))
    assert rc == 0 and " failures 0 " in out, out


def test_truncated_file(test_binary, tmp_path):
    path = tmp_path / "trunc.nc"
    _write(path, "NETCDF3_64BIT_DATA")
    size = path.stat().st_size
    with open(path, "r+b") as fh:
        fh.truncate(size - 30000)          # the last records' variables are cut
    # netCDF reads the cut variables as fill or fails; the test only needs cdf_read to refuse them,
    # so the reference read is skipped for refused ranges (see test_cdf_direct.F90)
    rc, out = _run(test_binary, path, _mode(test_binary, path))
    assert rc == 0 and " failures 0 " in out and " refused 0" not in out, out


def test_netcdf4_refused(test_binary, tmp_path):
    path = tmp_path / "t_nc4.nc"
    _write(path, "NETCDF4")
    rc, out = _run(test_binary, path, "direct")
    assert rc == 0 and "REFUSED" in out and "not a classic netCDF file" in out, out
