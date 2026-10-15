#!/usr/bin/env python3
"""
Script to convert ROMS `uscl` files to CDR forcing profiles.
Call with `-h` or `--help` for help.
"""

import argparse
import re
import numpy as np
import xarray as xr
from pathlib import Path

# Tracers always written to cdr_trcflx_profile (MARBL additionality vs ALT_CO2)
BASE_TRACER_NAMES = ("ALK", "DIC")

# `_uscl` variables holding the boundary fluxes of dedicated CDR_LITE tracers
CDR_ADD_VAR_RE = re.compile(
    r"^(?P<tracer>CDR_OAE_ALK\d+|CDR_OAE_DIC\d+|CDR_DOR_DIC\d+)"
    r"_add_(?P<bry>north|south|east|west)$"
)


def cdr_tracer_sort_key(name: str) -> tuple:
    """Sort key matching ROMS tracer order: OAE ALK/DIC pairs, then DOR DICs."""
    for group, pattern, slot in [(0, r"CDR_OAE_ALK(\d+)", 0),
                                 (0, r"CDR_OAE_DIC(\d+)", 1),
                                 (1, r"CDR_DOR_DIC(\d+)", 0)]:
        m = re.fullmatch(pattern, name)
        if m:
            return (group, int(m.group(1)), slot)
    return (2, 0, 0)


def discover_cdr_tracer_names(dataset: xr.Dataset, child_boundaries: dict) -> list:
    """CDR_LITE tracers whose fluxes are present on every active boundary."""
    active = {bry for bry, is_active in child_boundaries.items() if is_active}
    present = {}
    for varname in dataset.variables:
        m = CDR_ADD_VAR_RE.fullmatch(str(varname))
        if m:
            present.setdefault(m.group("tracer"), set()).add(m.group("bry"))
    return sorted((n for n, brys in present.items() if active and active <= brys),
                  key=cdr_tracer_sort_key)


class CDRUpscaler:
    """
    Class for transforming CDR-relevant boundary tracers to CDR forcing profiles.

    ALK and DIC additionality are always written to `cdr_trcflx_profile`. If the
    `_uscl` files also hold dedicated CDR_LITE tracers (`CDR_OAE_ALK*`,
    `CDR_OAE_DIC*`, `CDR_DOR_DIC*`), their fluxes follow along the tracer axis,
    and `cdr_trc_name` records which tracer each column is. ROMS uses the names
    to apply each column to the tracer of that name, skipping tracers the parent
    run does not have (e.g. ALK and DIC without MARBL).

    Use
    ---
    1. Initialize with a list of _joined_ `_uscl` files, e.g.
      `cu = CDRUpscaler(["roms_uscl.19920211170000.nc","roms_uscl.19920211180000.nc"])`
    This will open an xarray multi-file dataset with these files.

    2. A CDR forcing dataset based on the properties of these files can be initialized with
    `cu.create_cdr_dataset()`

    3. The CDR forcing dataset can then be populated based on the `uscl` data using
    `cu.populate_cdr_dataset()`

    4. The CDR forcing dataset can then be saved with
    `cu.save()`
    """

    def __init__(self,files: list):
        """
        Generate a CDRUpscaler instance from a list of joined, time-evolving `_uscl` files.

        Sets simulation- and domain-related attributes eta_rho, xi_rho, s_rho and time.
        Determines how many boundary columns are present to be transformed into CDR profiles.
        """

        self.files = files

        self.validate()

        print(f"Opening {len(self.files)} `_uscl` files...")
        self.uscl_dataset = xr.open_mfdataset(
            self.files,
            concat_dim="time",
            combine="nested",
            data_vars="minimal",
            coords="minimal"
        )
        self.eta_rho: xr.DataArray = self.uscl_dataset.eta_rho
        self.xi_rho: xr.DataArray = self.uscl_dataset.xi_rho
        self.s_rho: xr.DataArray = self.uscl_dataset.s_rho
        self.time: xr.DataArray = self.uscl_dataset.time

        # Determine active boundaries and add up all columns to get no. profiles
        self.child_boundaries: dict = {}
        self.n_profiles: int = 0
        for bry in ["north","south","east","west"]:
            v = f"ALK_add_{bry}"
            if v in self.uscl_dataset.variables:
                self.child_boundaries[bry] = True
                rho_var = "eta_rho" if bry.endswith("st") else "xi_rho"
                self.n_profiles = self.n_profiles + self.uscl_dataset.sizes.get(rho_var)
        print(f"Found {self.n_profiles} boundary tracer profiles to convert to CDR forcing profiles")

        cdr_tracer_names = discover_cdr_tracer_names(self.uscl_dataset, self.child_boundaries)
        self.tracer_names: list = list(BASE_TRACER_NAMES) + cdr_tracer_names
        if cdr_tracer_names:
            print("Found CDR tracers in `_uscl` files: " + ", ".join(cdr_tracer_names))

        # Initialize CDR dataset
        self.cdr_dataset: xr.Dataset | None = None

    def validate(self) -> None:
        """
        Run validation checks on input files
        """

        for i,f in enumerate(self.files):
            ds = xr.open_dataset(f)
            if i==0:
                reference_sizes = {k:v for k,v in ds.sizes.items() if k!="time"}
                continue
            ds_sizes = {k:v for k,v in ds.sizes.items() if k!="time"}
            if ds_sizes != reference_sizes:
                raise ValueError("Dimensions of input files do not match. "+
                                 "HINT: use ncjoin to join files before calling uscl_to_cdr")


    @property
    def filename_prefix(self) -> str:
        """
        Determine the common prefix of input `_uscl` files.

        Used to create default output filename for CDRUpscaler.save()

        Returns
        -------
        str: the common filename prefix of the CDRUpscaler inputs

        Examples
        --------
        CDRUpscaler(["output/roms_run_uscl.20250102130000.nc",
                     "output/roms_run_uscl.20250103130000.nc"
                    ]).filename_prefix()
        >> "output/roms_run"
        """

        prefixes = [f.split("_uscl")[0] for f in self.files]
        if not all([prefix == prefixes[0] for prefix in prefixes]):
            raise ValueError("Filenames do not share a common prefix")
        return prefixes[0]

    def create_cdr_dataset(self) -> None:
        """
        Initialize the CDR forcing dataset based on boundary geometry.

        Updates the CDRUpscaler.cdr_dataset attribute.

        Creates the following variables:
        - cdr_trcflx_profile: Series of time-evolving tracer flux profiles flux(time,depth,tracer,profile no.)
          ALK and DIC are the first two tracers; CDR_LITE tracers, if any, follow.
        - cdr_trc_name: Tracer names along the tracer axis (only with CDR_LITE tracers)
        - cdr_lon, cdr_lat: 1D arrays specifying where each tracer flux profile is located
        - cdr_layer_thickness: Series of time-evolving layer heights used to remap the profiles onto the parent grid
        """

        ds = xr.Dataset()
        n_trc = len(self.tracer_names)
        trc_dim = "two" if n_trc == 2 else "ntrc"

        ds["cdr_trcflx_profile"] = xr.DataArray(
            np.zeros((len(self.time), len(self.s_rho), n_trc, self.n_profiles)),
            dims=("cdr_time", "s_rho", trc_dim, "ncdr_prof"),
            attrs={
                "long_name": "tracer flux [mmol/s]",
                "units": "mmol/s",
            },
        )
        if n_trc > 2:
            ds["cdr_trc_name"] = xr.DataArray(
                np.array(self.tracer_names),
                dims=(trc_dim,),
                attrs={"long_name": "names of tracers stored in cdr_trcflx_profile"},
            )

        ds["cdr_time"] = xr.DataArray(
            self.uscl_dataset.ocean_time.values/86400.,
            dims=("cdr_time",),
            attrs={
                "long_name": "CDR forcing time",
                "units": "days since 2000/01/01",
            },
        )

        ds["cdr_lon"] = xr.DataArray(
            np.zeros(self.n_profiles),
            dims=("ncdr_prof",),
            attrs={
                "long_name": "longitude of CDR release [degrees East]",
                "units": "deg E",
            },
        )

        ds["cdr_lat"] = xr.DataArray(
            np.zeros(self.n_profiles),
            dims=("ncdr_prof",),
            attrs={
                "long_name": "latitude of CDR release [degrees North]",
                "units": "deg N",
            },
        )

        ds["cdr_layer_thickness"] = xr.DataArray(
            np.zeros((len(self.time), len(self.s_rho), self.n_profiles)),
            dims=("cdr_time", "s_rho", "ncdr_prof"),
            attrs={
                "long_name": "layer thicknesses of CDR release given in a vertical profile [m]",
                "units": "m",
            },
        )
        self.cdr_dataset = ds

    def populate_cdr_dataset(self):
        """
        Populates the CDRUpscaler.cdr_forcing() xarray Dataset created by CDRUpscaler.create_cdr_dataset().

        Remaps the time-evolving tracer flux columns of each boundary in the `_uscl` dataset to a depth
        profile in the CDR forcing dataset with a corresponding lat and lon.
        """


        def boundaries_to_profiles(var_prefix) -> np.ndarray:
            """Convert boundary columns to a sequence of profiles.

            Concatenates all boundary columns into a single list along
            the final axis (local lateral dimension).

            Note
            ----
            This formulation uses `.values` which triggers compute.
            See boundaries_to_profiles_xr for lazy implementation.
            """

            return np.concatenate(
                [self.uscl_dataset[f"{var_prefix}_{bry}"].values
                 for bry in ["north","east","south","west"]
                 if self.child_boundaries.get(bry)],
                axis=-1
            )

        def boundaries_to_profiles_xr(var_prefix) -> xr.DataArray:
            """boundaries_to_profiles using xarray/lazy concatenation.

            Should be faster, but isn't. Huge bottleneck during save()."""

            data_arrays=[]
            for bry in ["north","east","south","west"]:
                if self.child_boundaries.get(bry):
                    rho_var = "eta_rho" if bry.endswith("st") else "xi_rho"
                    da = self.uscl_dataset[f"{var_prefix}_{bry}"].rename({rho_var: "ncdr_prof"})
                    if "time" in da.dims:
                        da = da.rename({"time" : "cdr_time"})
                    data_arrays.append(da)
            return xr.concat(data_arrays, dim="ncdr_prof")

        print("Populating `cdr_lat`...")
        self.cdr_dataset["cdr_lat"][:] = boundaries_to_profiles("lat")
        print("Populating `cdr_lon`...")
        self.cdr_dataset["cdr_lon"][:] = boundaries_to_profiles("lon")
        for i_trc, name in enumerate(self.tracer_names):
            print(f"Populating `cdr_trcflx_profile` for tracer {name}...")
            self.cdr_dataset["cdr_trcflx_profile"][:,:,i_trc,:] = boundaries_to_profiles(f"{name}_add")
        print("Populating `cdr_layer_thickness`...")
        self.cdr_dataset["cdr_layer_thickness"][:] = boundaries_to_profiles("h")


    def save(self,filename: str | None = None) -> None:
        """
        Save the CDRForcing dataset to a netCDF file to be used in ROMS.

        Parameters
        ----------
        filename (str, optional): the filename of the saved netCDF file.
            Defaults to `CDR_forcing_from_<prefix>` where prefix is the
            shared prefix of the input file list found by
            `CDRUpscaler.filename_prefix()`
        """

        if not self.cdr_dataset:
            raise ValueError("No CDR forcing dataset to save. "+
                             "Call CDRUpscaler.create_cdr_dataset() and "+
                             "CDRUpscaler.populate_cdr_dataset() first")

        if not filename:
            basename = Path(self.filename_prefix).name
            parent = Path(self.filename_prefix).parent
            filepath = parent/f"cdr_release_profiles_from_{basename}.nc"
        else:
            filepath=Path(filename)

        print(f"Saving output to {filepath}")
        # ROMS (Fortran) reads tracer names as fixed-length characters
        # (NC_CHAR); xarray would otherwise write variable-length NC_STRING.
        encoding = {}
        if "cdr_trc_name" in self.cdr_dataset:
            encoding["cdr_trc_name"] = {"dtype": "S1", "char_dim_name": "trc_name_len"}
        self.cdr_dataset.to_netcdf(filepath, encoding=encoding)
        print("Done.")


def init_argparse() -> "ArgumentParser":
    """
    Initialize the ArgumentParser instance associated with this script.
    """

    description = """
    Utility for converting ROMS `<prefix>_uscl.<time>.nc` files to CDR forcing files.

    This program takes a series of `uscl` files from a nested/child simulation and
    creates a series of CDR forcing depth profiles that can be supplied to the parent
    simulation to propagate CDR-related tracers from child to parent conservatively.

    Each column at the edge of the child boundary is converted to a time-varying
    depth profile with lat/lon information such that ROMS can vertically remap it
    at the correct point on the parent grid.
    """

    parser = argparse.ArgumentParser(
        prog="uscl_to_cdr",
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("files", nargs="+", help="Input files")
    parser.add_argument("-o","--output",help="Output filename")

    return parser

if __name__ == "__main__":
    parser = init_argparse()
    args = parser.parse_args()

    cdr_upscaler = CDRUpscaler(sorted(args.files))
    cdr_upscaler.create_cdr_dataset()
    cdr_upscaler.populate_cdr_dataset()
    cdr_upscaler.save(args.output)
