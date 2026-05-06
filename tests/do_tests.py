import copy
import hashlib
import shutil
import json
import subprocess
import xarray as xr
from create_input_files import create_roms_inputs
from pathlib import Path
import f90nml as nl
import os

roms_root=Path(os.getenv("ROMS_ROOT"))

this_test = {}

#--------------------------------------------------------------------------------
# CPP key groups:
universal_cpp_keys = ["SOLVE3D","UV_ADV","UV_COR","UV_VIS2","TS_DIF2","MASKING","EXACT_RESTARTS"]
ocean_physics_cpp_keys = ["SALINITY","NONLIN_EOS","SPLIT_EOS",
                          "LMD_MIXING","LMD_KPP","LMD_BKPP","LMD_RIMIX","LMD_NONLOCAL"]
realistic_cpp_keys = ["CURVGRID","SPHERICAL","OBC_WEST","OBC_EAST","OBC_NORTH","OBC_SOUTH",
                      "OBC_M2FLATHER","OBC_M3ORLANSKI","OBC_TORLANSKI","Z_FRC_BRY",
                      "M2_FRC_BRY","M3_FRC_BRY","T_FRC_BRY"]
analytical_cpp_keys = ["ANA_GRID","ANA_INITIAL","ANA_SMFLUX","ANA_SRFLUX","ANA_STFLUX","ANA_SSFLUX"]
bgc_realistic_cpp_keys = ["BULK_FRC","SPONGE","SPONGE_TUNE","PCO2AIR_FORCING"]
marbl_cpp_keys = ["MARBL","MARBL_DIAGS","NOX_FORCING","NHY_FORCING"]
#--------------------------------------------------------------------------------

roms_nml=nl.read(roms_root/"src/namelist.nml")

def create_test_namelist_dict(input_dir):

    this_nml = copy.deepcopy(roms_nml)
    this_nml["TIME_STEPPING"]["ntimes"] = 20
    this_nml["TIME_STEPPING"]["dt"] = 40

    this_nml["GRID_SETTINGS"]["grdname"] = str(input_dir/"example_input_grid.nc")

    this_nml["S_COORD"]["theta_s"] = 6.0
    this_nml["S_COORD"]["theta_b"] = 6.0
    this_nml["S_COORD"]["hc"] = 25.0

    this_nml["PARAM_SETTINGS"]["NP_XI"] = 2
    this_nml["PARAM_SETTINGS"]["NP_ETA"] = 2
    this_nml["PARAM_SETTINGS"]["LLm"] = 39
    this_nml["PARAM_SETTINGS"]["MMm"] = 19
    this_nml["PARAM_SETTINGS"]["N"] = 10

    this_nml["INITIAL_CONDITIONS"]["ininame"] = str(input_dir/"example_input_bgc_initial_conditions.nc")

    this_nml["FORCING_FILES"]["frcfile"] = [
        str(input_dir/"example_input_boundary_forcing.nc"),
        str(input_dir/"example_input_surface_flux_forcing.nc"),
    ]

    this_nml["CALC_PFLX_SETTINGS"]["calc_pflx"] = False
    this_nml["CALC_PFLX_SETTINGS"]["timescale"] = 0

    this_nml["BASIC_OUTPUT_SETTINGS"]["wrt_file_his"] = True
    this_nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 800

    this_nml["BGC_SETTINGS"]["wrt_bgc_his"] = True
    this_nml["BGC_SETTINGS"]["interp_bgc_frc"] = True
    this_nml["BGC_SETTINGS"]["output_period_his"] = 400

    return this_nml

class ROMSConfiguration():
    def __init__(self,
                 cpp_keys: list[str],
                 location: Path
                 ):

        self.location = location
        self.cpp_keys = cpp_keys

        self.location.mkdir(parents=True)
        (self.location/"Makefile").symlink_to(roms_root/"Work/Makefile")

        with open(self.location/"cppdefs.opt","w") as f:
            for k in self.cpp_keys:
                f.write(f"#define {k}\n")
            f.write('#include "set_global_definitions.h"')

    def compile(self):
        subprocess.run(
            f"cd {self.location} && make BUILD_MODE=test",
            shell=True,
            text=True
        )

    def run(self, namelist_dict: dict):
        nl.write(
            namelist_dict,
            self.location/"namelist.nml",
            force=True)
        this_nml = nl.read(self.location/"namelist.nml")
        np_xi = this_nml["PARAM_SETTINGS"]["NP_XI"]
        np_eta = this_nml["PARAM_SETTINGS"]["NP_ETA"]
        ncpu = np_xi*np_eta

        subprocess.run(
            f"mpirun -n {ncpu} ./roms namelist.nml",
            cwd=self.location,
            shell=True,
            text=True
        )

def test_rivers_real(test_dir: Path, input_dir: Path) -> int:
    rivers_cpp = (universal_cpp_keys + ocean_physics_cpp_keys + realistic_cpp_keys +
                  ["SPONGE","SPONGE_TUNE"])

    rivers_conf = ROMSConfiguration(
        cpp_keys = rivers_cpp,
        location = test_dir
    )

    rivers_conf.compile()

    rivers_nml = create_test_namelist_dict(input_dir)
    rivers_nml["RIVER_FRC_SETTINGS"] = {
        "river_source": True,
        "river_analytical": False,
        "nriv": 1
    }
    rivers_nml["FORCING_FILES"]["frcfile"].append(str(input_dir/"example_input_river_forcing.nc"))
    rivers_conf.run(rivers_nml)

    # Check physics output
    pv = get_summary_value(test_dir, prefix="roms_his.20100101000000")
    this_test["rivers_real"]={"physics": pv}
    test_value=None
    with open(f"results_{environ}.json","r") as F:
        test_dict = json.load(F).get("rivers_real")
        if test_dict:
            test_value = test_dict.get("physics")

    return 0 if pv==test_value else 1

def test_rivers_ana(test_dir: Path, input_dir: Path) -> int:
    rivers_ana_cpp = (universal_cpp_keys + ocean_physics_cpp_keys + analytical_cpp_keys +
                      ["LMD_CONVEC","ANA_RIVER_USWC"])

    rivers_ana_conf = ROMSConfiguration(
        cpp_keys = rivers_ana_cpp,
        location = test_dir
    )
    additional_files_path=Path("additional_files").absolute()
    Path(test_dir/"ana_grid.h").symlink_to(additional_files_path/"rivers_ana_grid.h")
    Path(test_dir/"ana_init.h").symlink_to(additional_files_path/"rivers_ana_init.h")
    Path(test_dir/"ana_frc_river.h").symlink_to(additional_files_path/"rivers_ana_frc_river.h")

    rivers_ana_conf.compile()
    rivers_ana_nml = create_test_namelist_dict(input_dir)
    rivers_ana_nml["TIME_STEPPING"]["dt"] = 20
    rivers_ana_nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 400
    rivers_ana_nml["GRID_SETTINGS"] = {"grdname" : ""}
    rivers_ana_nml["INITIAL_CONDITIONS"] = {"ininame" : ""}
    rivers_ana_nml["FORCING_FILES"] = {"frcfile" : ""}
    rivers_ana_nml["PARAM_SETTINGS"].update({"NP_XI": 2, "NP_ETA": 2,"LLm" : 100, "MMm": 100, "N":10})
    rivers_ana_nml["RIVER_FRC_SETTINGS"] = {
        "river_source": True,
        "river_analytical": True,
        "nriv": 1
    }
    rivers_ana_conf.run(rivers_ana_nml)

    # Check output
    pv = get_summary_value(test_dir, prefix="roms_his.20000101000000")
    this_test["rivers_ana"] = {"physics": pv}
    test_value=None
    with open(f"results_{environ}.json","r") as F:
        test_dict = json.load(F).get("rivers_ana")
        if test_dict:
            test_value = test_dict.get("physics")
    return 0 if pv==test_value else 1

def test_pipes_real(test_dir: Path, input_dir: Path) -> int:
    pipes_cpp = (universal_cpp_keys + ocean_physics_cpp_keys + realistic_cpp_keys +
                 ["LMD_CONVEC"])

    pipes_conf = ROMSConfiguration(
        cpp_keys = pipes_cpp,
        location = test_dir
    )

    pipes_conf.compile()

    pipes_nml = create_test_namelist_dict(input_dir)
    pipes_nml["TIME_STEPPING"]["dt"] = 20
    pipes_nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 400
    pipes_nml["FORCING_FILES"]["frcfile"].append(str(input_dir/"example_input_pipe_forcing.nc"))
    pipes_nml["PIPE_FRC_SETTINGS"] = {
        "pipe_source": True,
        "p_analytical": False,
        "npip": 1
    }

    pipes_conf.run(pipes_nml)

    # Check physics output
    pv = get_summary_value(test_dir, prefix="roms_his.20100101000000")
    this_test["pipes_real"]={"physics": pv}
    test_value=None
    with open(f"results_{environ}.json","r") as F:
        test_dict = json.load(F).get("pipes_real")
        if test_dict:
            test_value = test_dict.get("physics")

    return 0 if pv==test_value else 1

def test_pipes_ana(test_dir: Path, input_dir: Path) -> int:
    pipes_ana_cpp = (universal_cpp_keys + ocean_physics_cpp_keys + analytical_cpp_keys +
                      ["LMD_CONVEC","ANA_PIPES"])

    pipes_ana_conf = ROMSConfiguration(
        cpp_keys = pipes_ana_cpp,
        location = test_dir
    )
    additional_files_path=Path("additional_files").absolute()
    Path(test_dir/"ana_grid.h").symlink_to(additional_files_path/"pipes_ana_grid.h")
    Path(test_dir/"ana_init.h").symlink_to(additional_files_path/"pipes_ana_init.h")
    Path(test_dir/"ana_pipe_frc.h").symlink_to(additional_files_path/"pipes_ana_pipe_frc.h")

    pipes_ana_conf.compile()
    pipes_ana_nml = create_test_namelist_dict(input_dir)
    pipes_ana_nml["TIME_STEPPING"]["dt"] = 60
    pipes_ana_nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 1200
    pipes_ana_nml["GRID_SETTINGS"] = {"grdname" : ""}
    pipes_ana_nml["INITIAL_CONDITIONS"] = {"ininame" : ""}
    pipes_ana_nml["FORCING_FILES"] = {"frcfile" : ""}
    pipes_ana_nml["PARAM_SETTINGS"].update({"NP_XI":2, "NP_ETA":2, "LLm" : 100, "MMm": 100, "N":10})
    pipes_ana_nml["PIPE_FRC_SETTINGS"] = {
        "pipe_source": True,
        "p_analytical": True,
        "npip": 1
    }
    pipes_ana_conf.run(pipes_ana_nml)

    # Check output
    pv = get_summary_value(test_dir, prefix="roms_his.20000101000000")
    this_test["pipes_ana"] = {"physics": pv}
    test_value=None
    with open(f"results_{environ}.json","r") as F:
        test_dict = json.load(F).get("pipes_ana")
        if test_dict:
            test_value = test_dict.get("physics")
    return 0 if pv==test_value else 1

def test_flux_frc(test_dir: Path, input_dir: Path) -> int:
    flux_frc_cpp = (universal_cpp_keys + ocean_physics_cpp_keys + realistic_cpp_keys +
                    ["LMD_CONVEC","ADV_ISONEUTRAL"])

    flux_frc_conf = ROMSConfiguration(
        cpp_keys = flux_frc_cpp,
        location = test_dir
    )

    flux_frc_conf.compile()
    flux_frc_nml = create_test_namelist_dict(input_dir)
    flux_frc_nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 400
    flux_frc_conf.run(flux_frc_nml)

    pv = get_summary_value(test_dir, prefix="roms_his.20100101000000")
    this_test["flux_frc"] = {"physics": pv}
    test_value=None
    with open(f"results_{environ}.json","r") as F:
        test_dict = json.load(F).get("flux_frc")
        if test_dict:
            test_value = test_dict.get("physics")
    return 0 if pv==test_value else 1

def test_filament(test_dir: Path, input_dir: Path) -> int:
    filament_cpps = (universal_cpp_keys + analytical_cpp_keys +
                     ["ANA_VMIX","NS_PERIODIC","EW_PERIODIC","FILAMENT_IDEAL"])

    filament_conf = ROMSConfiguration(
        cpp_keys = filament_cpps,
        location = test_dir
    )
    additional_files_path=Path("additional_files").absolute()
    Path(test_dir/"ana_grid.h").symlink_to(additional_files_path/"filament_grid.h")
    Path(test_dir/"ana_init.h").symlink_to(additional_files_path/"filament_init.h")
    Path(test_dir/"ana_vmix.h").symlink_to(additional_files_path/"filament_vmix.h")


    filament_conf.compile()
    filament_nml = create_test_namelist_dict(input_dir)
    filament_nml["TIME_STEPPING"].update({"dt": 5, "ndtfast": 60})
    filament_nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 50
    filament_nml["GRID_SETTINGS"] = {"grdname" : ""}
    filament_nml["INITIAL_CONDITIONS"] = {"ininame" : ""}
    filament_nml["FORCING_FILES"] = {"frcfile" : ""}
    filament_nml["PARAM_SETTINGS"].update({"NP_XI":2,"NP_ETA":2,"LLm" : 64, "MMm": 64, "N":32})
    filament_nml["VERTICAL_MIXING_SETTINGS"] = {"akv_bak" : 0., "akt_bak": 0.}
    filament_nml["LIN_RHO_EOS_SETTINGS"] = {"Tcoef" : 0.2, "T0": 1.0}
    filament_nml["RHO0_SETTINGS"]["rho0"] = 1000
    filament_nml["S_COORD"]["theta_b"] = 2.0
    filament_conf.run(filament_nml)

    pv = get_summary_value(test_dir, prefix="roms_his.20000101000000")
    this_test["filament"] = {"physics": pv}
    test_value=None
    with open(f"results_{environ}.json","r") as F:
        test_dict = json.load(F).get("filament")
        if test_dict:
            test_value = test_dict.get("physics")
    return 0 if pv==test_value else 1

def test_cdr(test_dir: Path, input_dir: Path) -> int:
    cdr_cpps = (universal_cpp_keys + ocean_physics_cpp_keys +
                realistic_cpp_keys + bgc_realistic_cpp_keys +
                marbl_cpp_keys + ["CDR_FORCING"])

    cdr_conf = ROMSConfiguration(
        cpp_keys = cdr_cpps,
        location = test_dir
    )

    cdr_conf.compile()

    cdr_nml = create_test_namelist_dict(input_dir)
    cdr_nml["TIME_STEPPING"]["dt"] = 40
    cdr_nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 400
    cdr_nml["BGC_SETTINGS"]["output_period_his"] = 400
    cdr_nml["FORCING_FILES"] = {
        "frcfile": [
            str(input_dir/"example_input_boundary_forcing.nc"),
            str(input_dir/"example_input_surface_forcing.nc"),
            str(input_dir/"example_input_bgc_surface_forcing.nc"),
            str(input_dir/"example_input_river_forcing.nc"),
            str(input_dir/"example_input_bgc_boundary_forcing.nc")
        ]}
    cdr_nml["PARAM_SETTINGS"]["ntrc_bio"] = 32

    counter = 0

    # CDR_3d test
    subprocess.run("rm *.nc",cwd=test_dir,shell=True,text=True)
    cdr_3d_nml = copy.deepcopy(cdr_nml)
    cdr_3d_nml["FORCING_FILES"]["frcfile"] += [str(input_dir/"cdr_forcing_3d.nc")]
    cdr_3d_nml["CDR_FRC_SETTINGS"] = {"cdr_source": True,
                                      "forcing_3d": True,
                                      "relocate_to_wet_pts": True}
    cdr_conf.run(cdr_3d_nml)

    pv = get_summary_value(test_dir, prefix="roms_his.20100101000000")
    bv = get_summary_value(test_dir, prefix="roms_bgc.20100101000000")

    this_test["cdr_3d"] = {"physics": pv, "bgc": bv}
    test_value_phys=None
    test_value_bgc=None

    with open(f"results_{environ}.json","r") as F:
        test_dict = json.load(F).get("cdr_3d")
        if test_dict:
            test_value_phys = test_dict.get("physics")
            test_value_bgc = test_dict.get("bgc")

    retval = 0 if ( (test_value_phys == pv) and (test_value_bgc == bv)) else 1
    counter = counter+retval

    # CDR_dp test
    subprocess.run("rm *.nc",cwd=test_dir,shell=True,text=True)
    cdr_dp_nml = copy.deepcopy(cdr_nml)
    cdr_dp_nml["BGC_SETTINGS"]["output_period_his"] = 200
    cdr_dp_nml["CDR_FRC_SETTINGS"].update({"cdr_source": True,
                                           "forcing_depth_profiles": True,
                                           "cdr_file": str(input_dir/"cdr_forcing_dp.nc"),
                                           "relocate_to_wet_pts": True,
                                           "nz_chd" : 10})
    cdr_conf.run(cdr_dp_nml)

    pv = get_summary_value(test_dir, prefix="roms_his.20100101000000")
    bv = get_summary_value(test_dir, prefix="roms_bgc.20100101000000")

    this_test["cdr_dp"] = {"physics": pv, "bgc": bv}
    test_value_phys=None
    test_value_bgc=None

    with open(f"results_{environ}.json","r") as F:
        test_dict = json.load(F).get("cdr_dp")
        if test_dict:
            test_value_phys = test_dict.get("physics")
            test_value_bgc = test_dict.get("bgc")

    retval = 0 if ( (test_value_phys == pv) and (test_value_bgc == bv)) else 1
    counter = counter+retval

    # # CDR_parm test
    subprocess.run("rm *.nc",cwd=test_dir,shell=True,text=True)
    cdr_parm_nml = copy.deepcopy(cdr_nml)
    cdr_parm_nml["CDR_FRC_SETTINGS"].update({"cdr_source": True,
                                             "forcing_parameterized": True,
                                             "relocate_to_wet_pts": True,
                                             "cdr_file": str(input_dir/"cdr_forcing_parm.nc")})
    cdr_conf.run(cdr_parm_nml)
    pv = get_summary_value(test_dir, prefix="roms_his.20100101000000")
    bv = get_summary_value(test_dir, prefix="roms_bgc.20100101000000")

    this_test["cdr_parameterized"] = {"physics": pv, "bgc": bv}
    test_value_phys=None
    test_value_bgc=None

    with open(f"results_{environ}.json","r") as F:
        test_dict = json.load(F).get("cdr_parameterized")
        if test_dict:
            test_value_phys = test_dict.get("physics")
            test_value_bgc = test_dict.get("bgc")

    retval = 0 if ( (test_value_phys == pv) and (test_value_bgc == bv)) else 1
    counter = counter+retval

    return counter

def test_bgc_real_bec(test_dir: Path, input_dir: Path) -> int:
    bgc_real_BEC_cpps = (universal_cpp_keys + ocean_physics_cpp_keys +
                          realistic_cpp_keys + bgc_realistic_cpp_keys +
                          ["TIDES","BIOLOGY_BEC2"])

    bec_conf = ROMSConfiguration(
        cpp_keys = bgc_real_BEC_cpps,
        location = test_dir
    )

    bec_conf.compile()
    bec_nml = create_test_namelist_dict(input_dir)
    bec_nml["TIME_STEPPING"]["dt"] = 20
    bec_nml["TIME_STEPPING"]["ntimes"] = 10
    bec_nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 100
    bec_nml["BGC_SETTINGS"]["output_period_his"] = 400
    bec_nml["PARAM_SETTINGS"]["ntrc_bio"] = 26
    bec_nml["TIDES_SETTINGS"].update({"pot_tides": True, "ntides": 2})
    bec_nml["FORCING_FILES"] = {
        "frcfile": [
            str(input_dir/"example_input_boundary_forcing.nc"),
            str(input_dir/"example_input_surface_forcing.nc"),
            str(input_dir/"example_input_bgc_surface_forcing.nc"),
            str(input_dir/"example_input_river_forcing.nc"),
            str(input_dir/"example_input_bgc_boundary_forcing.nc"),
            str(input_dir/"example_input_tides.nc")
        ]}

    bec_conf.run(bec_nml)

    pv = get_summary_value(test_dir, prefix="roms_his.20100101000000")
    bv = get_summary_value(test_dir, prefix="roms_bgc.20100101000000")

    this_test["bgc_real_bec"] = {"physics": pv, "bgc": bv}
    test_value_phys=None
    test_value_bgc=None

    with open(f"results_{environ}.json","r") as F:
        test_dict = json.load(F).get("bgc_real_bec")
        if test_dict:
            test_value_phys = test_dict.get("physics")
            test_value_bgc = test_dict.get("bgc")

    return 0 if ( (test_value_phys == pv) and (test_value_bgc == bv)) else 1

def test_bgc_real_marbl(test_dir: Path, input_dir: Path) -> int:
    bgc_real_marbl_cpps = ( universal_cpp_keys + ocean_physics_cpp_keys +
                            realistic_cpp_keys + bgc_realistic_cpp_keys +
                            marbl_cpp_keys + ["TIDES"])

    marbl_conf = ROMSConfiguration(
        cpp_keys = bgc_real_marbl_cpps,
        location = test_dir
    )

    marbl_conf.compile()
    additional_files_path=Path("additional_files").absolute()
    Path(test_dir/"marbl_in").symlink_to(additional_files_path/"marbl_in")
    marbl_nml = create_test_namelist_dict(input_dir)
    marbl_nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 100
    marbl_nml["BGC_SETTINGS"]["output_period_his"] = 100
    marbl_nml["MARBL_BIOGEOCHEMISTRY_SETTINGS"]["marbl_diagnostics_to_write"] = ["DCO2STAR_ALT_CO2","graze_diat_zint","FESEDFLUX"]
    marbl_nml["TIME_STEPPING"].update({"dt":  20, "ntimes": 10})
    marbl_nml["TIDES_SETTINGS"].update({"pot_tides": True, "ntides": 2})
    marbl_nml["PARAM_SETTINGS"]["ntrc_bio"] = 32
    marbl_nml["FORCING_FILES"] = {
        "frcfile": [
            str(input_dir/"example_input_boundary_forcing.nc"),
            str(input_dir/"example_input_surface_forcing.nc"),
            str(input_dir/"example_input_bgc_surface_forcing.nc"),
            str(input_dir/"example_input_river_forcing.nc"),
            str(input_dir/"example_input_bgc_boundary_forcing.nc"),
            str(input_dir/"example_input_tides.nc")
        ]}

    marbl_conf.run(marbl_nml)

    pv = get_summary_value(test_dir, prefix="roms_his.20100101000000")
    bv = get_summary_value(test_dir, prefix="roms_bgc.20100101000000")

    this_test["bgc_real_marbl"] = {"physics": pv, "bgc": bv}
    test_value_phys=None
    test_value_bgc=None

    with open(f"results_{environ}.json","r") as F:
        test_dict = json.load(F).get("bgc_real_marbl")
        if test_dict:
            test_value_phys = test_dict.get("physics")
            test_value_bgc = test_dict.get("bgc")

    return 0 if ( (test_value_phys == pv) and (test_value_bgc == bv)) else 1


def get_summary_value(test_dir: Path, prefix: str):
    test_dir = test_dir.absolute()
    subprocess.run(
        f"ncjoin {prefix}.?.nc",
        cwd=test_dir,
        text=True, shell=True
    )
    val = 0
    ds=xr.open_dataset(f"{test_dir}/{prefix}.nc")

    h = hashlib.sha256()
    for name in sorted(ds.data_vars):  # sort for determinism
        h.update(ds[name].values.tobytes())
    return h.hexdigest()



if __name__ == "__main__":
    environ = "laptop"
    my_dir = Path("tmp")
    inp_dir = Path("tmp2")
    inp_dir.mkdir()
    inp_dir=inp_dir.absolute()
    create_roms_inputs(inp_dir)
    counter = 0
    failures = []
    # Rivers Real
    result  = test_rivers_real(
        test_dir = my_dir/"rivers_real",
        input_dir = inp_dir
    )
    if result != 0: failures.append("rivers_real")
    # Rivers Ana
    result = test_rivers_ana(
        test_dir = my_dir/"rivers_ana",
        input_dir = inp_dir
    )
    if result != 0: failures.append("rivers_ana")
    # Pipes Real
    result = test_pipes_real(
        test_dir = my_dir/"pipes_real",
        input_dir = inp_dir
    )
    if result != 0: failures.append("pipes_real")
    # Pipes Ana
    result = test_pipes_ana(
        test_dir = my_dir/"pipes_ana",
        input_dir = inp_dir
    )
    if result != 0: failures.append("pipes_ana")
    # Flux frc
    result = test_flux_frc(
        test_dir = my_dir/"flux_frc",
        input_dir = inp_dir
    )
    if result != 0: failures.append("flux_frc")
    #Filament
    result = test_filament(
        test_dir = my_dir/"filament",
        input_dir = inp_dir
    )
    if result != 0: failures.append("filament")
    #CDR
    result = test_cdr(
        test_dir = my_dir/"cdr",
        input_dir = inp_dir
    )
    if result != 0: failures.append("cdr")
    # bgc real BEC
    result = test_bgc_real_bec(
        test_dir = my_dir/"bgc_real_bec",
        input_dir = inp_dir
    )
    if result != 0: failures.append("bgc_real_bec")
    # bgc real MARBL
    result = test_bgc_real_marbl(
        test_dir = my_dir/"bgc_real_marbl",
        input_dir = inp_dir
    )
    if result != 0: failures.append("bgc_real_marbl")

    if len(failures) != 0:
        print(f"THE FOLLOWING TESTS FAILED: {failures}")
        print(f"If you understand why and would like to update results_{environ}.json, paste:")
        print(json.dumps(this_test))
    else:
        print("ALL TESTS PASSED")

