#!/usr/bin/env python3
"""
ROMS simulation output log parser returning json-format run info.

Usage
-----
Call using `parse_log my_output_log.txt` to produce output `my_output_log_<identifier>.json`,
where <identifier> is a SLURM jobID or the current date.

Details
-------
This script provides a modular parser for ROMS stdout files.
Each section of the ROMS output is parsed by a different function,
decorated with the `@register_logparser("attr0",...)` decorator.
A class `ROMSSimulationLog` is populated by these functions.
Applying the decorator to a new function causes the script to pass the output
file through that function to extract any relevant info and append it to the class.
The decorator should specify any attrs on the class it intends to define.
The script then dumps all the collected info to a json file with a
unique identifier (either the SLURM jobID or current date.)
"""

import os
import re
import sys
import json
import inspect
import datetime
import subprocess
from pathlib import Path
from typing import Callable

# Initialize list of parsers registered by the `register_logparser` decorator:
_parser_registry: list[tuple[Callable, tuple[str, ...]]] = []

def register_logparser(*attribute_names: str):
    """Decorator to register a parser function and the attributes it sets"""
    def decorator(func: Callable) -> Callable:
        _parser_registry.append((func, attribute_names))
        return func
    return decorator

class ROMSSimulationLog:
    def __init__(self,
                 log_file: Path):
        """Initialize this ROMSSimulationLog instance.

        Population is in three steps:
        1. parses the provided log file.
        2. query SLURM for additional info using `sacct` on the jobID extracted from the file.
        3. Determine the calling machine's name from an environment variable.
        """

        self.log_file = log_file
        with open(self.log_file, "r") as f:
            self.lines = f.readlines()
        # Automatically initialize all attributes declared by parsers
        for _, attrs in _parser_registry:
            for attr in attrs:
                if not hasattr(self, attr):
                    if (
                            attr.endswith("s") or
                            attr.endswith("_files") or
                            attr.startswith("output_variables")
                    ):
                        setattr(self, attr, [])
                    else:
                        setattr(self, attr, None)

        self._parse_logfile()

        # Query most remaining attrs from Slurm
        self._query_job_id()

        # Lastly, get machine info
        rcac = os.environ.get("RCAC_CLUSTER", "").lower()
        lmod = os.environ.get("LMOD_SYSHOST", "").lower()
        if rcac == "anvil":
            self.machine = "anvil"
        elif lmod == "perlmutter":
            self.machine = "perlmutter"
        else:
            self.machine = None


    def _parse_logfile(self):
        """Read the logfile by line, passing each line to each parser.

        The parser returns a new line to move to.
        """
        current_line = 0
        while current_line < len(self.lines):
            for parser_func, _ in _parser_registry:
                next_line = parser_func(self, current_line)
                if next_line != current_line:
                    current_line = next_line
                    break
            else:
                current_line += 1

    @classmethod
    def list_attributes(cls) -> list[str]:
        parsed_attrs = [attr for _, attrs in _parser_registry for attr in attrs]
        # slurm_attrs = ["slurm_maxrss","slurm_elapsed","slurm_state","slurm_exitcode","slurm_totalcpu"]
        slurm_attrs = [
            "slurm_maxrss","slurm_elapsed","slurm_state","slurm_exitcode","slurm_totalcpu",
            "slurm_job_name","slurm_user","slurm_partition","slurm_start_time","slurm_end_time",
        ]
        return parsed_attrs + slurm_attrs + ["machine",]

    @property
    def ntracers(self):
        """The number of tracers in the simulation."""
        return len(self.tracers)

    @property
    def walltime(self):
        """The walltime determined within ROMS."""
        return self.timestep_walltime[-1]

    @property
    def npoints(self):
        """The number of gridpoints in the ROMS domain."""
        return self.nx * self.ny * self.nz

    @property
    def performance_number(self):
        """Diagnostic to estimate system performance based on other attrs."""
        return (self.ncpus * self.walltime) / (self.ntimes * self.npoints)

    def _query_job_id(self):
        """Determine further information from SLURM's sacct function.

        Uses the `job_id` attr extracted from `_parse_logfile`
        """
        if not self.job_id:
            return
        try:
            result = subprocess.run(
                [
                    "sacct", "-j", str(self.job_id),
                    "--format=JobIDRaw,JobName,User,Partition,MaxRSS,Elapsed,Start,End,State,ExitCode,TotalCPU,Nodelist",
                    "--parsable2", "--noheader"
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
        except subprocess.CalledProcessError as e:
            print(f"[slurm query failed] {e.stderr.strip()}")
            return

        lines = result.stdout.strip().splitlines()
        if not lines:
            return

        top_level = None
        batch_step = None

        for line in lines:
            parts = line.strip().split('|')
            if len(parts) < 6:
                continue  # malformed line

            jobid_raw = parts[0]

            if jobid_raw.endswith(".batch"):
                batch_step = parts
            elif '.' not in jobid_raw:
                top_level = parts

        job_info = batch_step or top_level

        # if we have both, merge them; otherwise fall back to whichever we have
        if top_level and batch_step:
            # overwrite user/jobname/partition from top-level
            job_info[1:4] = top_level[1:4]
        if not job_info:
            return  # nothing usable found

        (
            jobid_raw, job_name, user, partition,
            maxrss, elapsed, start, end, state, exitcode, totalcpu, nodelist
        ) = job_info[:12]

        self.slurm_job_name = job_name or None
        self.slurm_user = user or None
        self.slurm_partition = partition or None
        self.slurm_start_time = start or None
        self.slurm_end_time = end or None
        self.slurm_maxrss = maxrss or None
        self.slurm_elapsed = elapsed or None
        self.slurm_state = state or None
        self.slurm_exitcode = exitcode or None
        self.slurm_totalcpu = totalcpu or None


    def _to_serializable(self,obj):
        """recursively make sure Paths etc become plain serializable types"""
        if isinstance(obj, Path):
            return str(obj)
        if isinstance(obj, dict):
            return {k: self._to_serializable(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [self._to_serializable(v) for v in obj]
        return obj

    def to_dict(self):
        """Return all serializable attributes and properties."""
        attrs = {}

        # include normal attributes
        for name in self.list_attributes():
            if hasattr(self, name):
                attrs[name] = self._to_serializable(getattr(self, name))

        # include class-level @property attributes
        for name, attr in inspect.getmembers(type(self)):
            if isinstance(attr, property) and not name.startswith("_"):
                try:
                    attrs[name] = self._to_serializable(getattr(self, name))
                except Exception:
                    pass  # ignore properties that fail to compute

        return attrs

    def to_json(self, path: Path | None = None, indent: int = 2):
        """Dump the output of `to_dict` to json format, optionally to file."""
        data = json.dumps(self.to_dict(), indent=indent)
        if path:
            Path(path.parent).mkdir(parents=True, exist_ok=True)
            Path(path).write_text(data)
        return data



@register_logparser("git_hash")
def parse_git_hash(log, i):
    """Get the git hash of the ROMS repo from the logfile."""
    line = log.lines[i].strip().lower()
    if "git hash" in line:
        log.git_hash = line.split(":")[-1].strip()
        return i+1
    return i

@register_logparser("job_id","job_task_id")
def parse_job_id(log, i):
    """Get the SLURM or PBS job ID from the logfile."""
    line = log.lines[i].strip()

    # job id: SLURM Job ID: <id>   OR   PBS Job ID: <id>
    if line.startswith(("SLURM Job ID:", "PBS Job ID:")):
        log.job_id = line.split(":", 1)[1].strip()
        return i+1

    # array id: SLURM Array Task ID: <id>   OR   PBS Array Index: <id>
    if line.startswith(("SLURM Array Task ID:", "PBS Array Index:")):
        log.job_array_task_id = line.split(":", 1)[1].strip()
        return i+1

    return i

@register_logparser("cppdefs")
def parse_cppdefs(log: ROMSSimulationLog, current_line: int) -> int:
    """Get the list of C preprocessor keys from the logfile."""
    line = log.lines[current_line]
    if "<cppdefs.opt>" not in line:
        return current_line
    j = current_line + 1
    while j < len(log.lines):
        l = log.lines[j].strip()
        if not l or (l.startswith("<") and l.endswith(">")):
            break
        log.cppdefs.extend(l.split())
        j += 1
    return j


KV_PARAMS = {"ntimes", "dt", "ndtfast", "ninfo", "theta_s", "theta_b", "hc", "rho0","visc2"}

@register_logparser(*KV_PARAMS)
def parse_key_value_params(log: "ROMSSimulationLog", current_line: int) -> int:
    """Parse a selection of related parameters formatted as `key = value` in the logfile."""

    line = log.lines[current_line]
    if "=" not in line:
        return current_line

    parts = line.split("=", 1)
    key = parts[0].strip()
    if key not in KV_PARAMS:
        return current_line

    rhs = parts[1].strip()
    token = rhs.split()[0]

    try:
        val = float(token) if "." in token or "E" in token.upper() else int(token)
    except ValueError:
        val = token

    setattr(log, key, val)
    return current_line + 1


@register_logparser("grid_file")
def parse_grid_file(log: "ROMSSimulationLog", current_line: int) -> int:
    """Parse the grid file path from the logfile."""
    line = log.lines[current_line]
    if not line.strip().startswith("grid file:"):
        return current_line
    log.grid_file = Path(line.split(":", 1)[1].strip())
    return current_line + 1


@register_logparser("forcing_files")
def parse_forcing_files(log: "ROMSSimulationLog", current_line: int) -> int:
    """Parse the list of forcing file paths from the logfile."""
    line = log.lines[current_line]
    if not line.strip().startswith("forcing data file(s):"):
        return current_line

    files = [Path(line.split(":", 1)[1].strip())]
    j = current_line + 1
    while j < len(log.lines):
        nxt = log.lines[j]
        stripped = nxt.strip()
        if not stripped or not nxt.startswith(" " * 2):
            break
        files.append(Path(stripped))
        j += 1
    log.forcing_files = files
    return j


@register_logparser("initial_condition_rec","initial_condition_file")
def parse_initial_condition(log: "ROMSSimulationLog", current_line: int) -> int:
    """Parse the initial condition file path from the logfile."""
    line = log.lines[current_line]
    if not line.strip().startswith("initial condition"):
        return current_line

    tokens = line.split()
    if "rec" in tokens:
        idx = tokens.index("rec")
        log.initial_condition_rec = int(tokens[idx + 2])  # skip '=' then value
    if "file" in tokens:
        idx = tokens.index("file")
        log.initial_condition_file = Path(tokens[idx + 2].strip("'"))
    return current_line + 1

@register_logparser("ncpus","np_xi","np_eta","nx","ny","nz")
def parse_node_grid_info(log: "ROMSSimulationLog", current_line: int) -> int:
    """Parse the CPU distribution info from the logfile."""
    line = log.lines[current_line].strip()

    if not line.startswith("NUMBER OF NODES:"):
        return current_line

    # compile locally (readable, one and done)
    nodes_re = re.compile(r"NUMBER OF NODES:\s*(?P<ncpus>\d+)")
    procs_re = re.compile(r"\((?P<np_xi>\d+)\s*x\s*(?P<np_eta>\d+)\)")
    grid_re = re.compile(r"GRID:\s*(?P<nx>\d+)\s*x\s*(?P<ny>\d+)\s*x\s*(?P<nz>\d+)")

    m = nodes_re.search(line)
    if m:
        log.ncpus = int(m.group("ncpus"))

    m = procs_re.search(line)
    if m:
        log.np_xi  = int(m.group("np_xi"))
        log.np_eta = int(m.group("np_eta"))

    m = grid_re.search(line)
    if m:
        log.nx = int(m.group("nx"))
        log.ny = int(m.group("ny"))
        log.nz = int(m.group("nz"))

    return current_line + 1

@register_logparser("tracers")
def parse_tracers(log: "ROMSSimulationLog", current_line: int) -> int:
    """Parse the list of tracers from the logfile."""
    line = log.lines[current_line].strip()
    if not line.startswith("TRACER NO.:"):
        return current_line

    if not hasattr(log, "tracers"):
        log.tracers = []

    j = current_line
    while j < len(log.lines):
        line = log.lines[j].strip()
        if not line.startswith("TRACER NO.:"):
            break

        # SHORT NAME is always the next line
        short_line = log.lines[j + 1].strip()
        if short_line.startswith("SHORT NAME:"):
            short_name = short_line.split(":", 1)[1].strip()
            log.tracers.append(short_name)

        # skip ahead until the end of this block (the dashed line)
        while j < len(log.lines) and not log.lines[j].strip().startswith("-----------"):
            j += 1
        j += 1  # move past the dashes

    return j

################################################################################
# OUTPUT SETTINGS:
@register_logparser("output_period_rst","nrpf_rst")
def parse_rst(log, i):
    """ Parse the restart file output settings from the logfile."""
    if "ocean_vars :: restart file" not in log.lines[i]:
        return i

    line = log.lines[i]

    # recs/file
    m = re.search(r"recs/file\s*=\s*(\d+)", line)
    log.nrpf_rst = int(m.group(1)) if m else None

    # period or monthly
    if "monthly_restarts" in line:
        log.output_period_rst = "monthly"
    else:
        m = re.search(r"output_period_rst\s*=\s*([0-9.]+)", line)
        log.output_period_rst = float(m.group(1)) if m else None

    return i + 1


# Generic parsing function for bgc/phys/cstar:
def _parse_output_block(log, start_i: int,
                        period_attr: str,
                        nrpf_attr: str,
                        vars_attr: str | None):
    """
    Helper function to parse related output settings sections from the logfile.

    Determines shared values:
    - period_attr: attribute on log to store output period (float or "monthly")
    - nrpf_attr: attribute on log to store recs/file (int)
    - vars_attr: attribute to store list of var names (list[str]) — or None to skip
    returns new index
    """
    line0 = log.lines[start_i]

    # detect period, nrpf, monthly in header line
    monthly = False
    period = None
    nrpf = None

    # recs/file
    m = re.search(r"recs/file\s*=\s*(\d+)", line0)
    if m: nrpf = int(m.group(1))

    # period
    m = re.search(r"output_period\s*=\s*([0-9.]+)", line0)
    if m: period = float(m.group(1))

    # monthly flags
    if "monthly_averages" in line0:
        monthly = True

    # assign period/nrpf on the spot
    setattr(log, period_attr, "monthly" if monthly else period)
    setattr(log, nrpf_attr, nrpf)

    # scan down to first dashed line
    dash_re = re.compile(r"^-{4,}$")  # 4+ hyphens
    i = start_i + 1
    while i < len(log.lines) and not dash_re.match(log.lines[i].strip()):
        i += 1
    # skip 2+ dash header block lines
    seen = 0
    while i < len(log.lines) and dash_re.match(log.lines[i].strip()) and seen < 2:
        seen += 1
        i += 1
    # skip optional header titles
    while i < len(log.lines) and not dash_re.match(log.lines[i].strip()):
        i += 1
    # this dashed row ends header
    while i < len(log.lines) and dash_re.match(log.lines[i].strip()):
        i += 1

    # now variable rows until next dashed
    vars = []
    while i < len(log.lines):
        s = log.lines[i].strip()
        if dash_re.match(s):
            break
        # if row looks like "<name> ... T ..."
        if vars_attr:
            parts = s.split()
            if len(parts) >= 2 and parts[1].upper().startswith("T"):
                vars.append(parts[0])
        i += 1

    if vars_attr and vars:
        setattr(log, vars_attr, vars)

    # skip final dashed lines
    while i < len(log.lines) and dash_re.match(log.lines[i].strip()):
        i += 1

    return i


# ------------------- ocean vars -------------------

@register_logparser("output_period_phys_his","nrpf_phys_his","output_variables_phys")
def parse_phys_his(log, i):
    """Parse the 'history' file (instantaneous output) settings from the logfile.
    """
    if "ocean_vars :: history file" not in log.lines[i]:
        return i
    return _parse_output_block(
        log, i,
        period_attr="output_period_phys_his",
        nrpf_attr="nrpf_phys_his",
        vars_attr="output_variables_phys"
    )

@register_logparser("output_period_phys_avg","output_nrpf_phys_avg","output_variables_phys")
def parse_phys_avg(log, i):
    """Parse the averaged output file settings from the logfile."""
    if "ocean_vars :: average file" not in log.lines[i]:
        return i
    return _parse_output_block(
        log, i,
        period_attr="output_period_phys_avg",
        nrpf_attr="output_nrpf_phys_avg",
        vars_attr="output_variables_phys"
    )

# ------------------- bgc -------------------

@register_logparser("output_period_bgc_his", "output_nrpf_bgc_his", "output_variables_bgc")
def parse_bgc_his(log, i):
    """Parse the biogeochemistry history file (instantaneous output) settings from the logfile."""
    if "bgc :: history file" not in log.lines[i].strip().lower():
        return i
    return _parse_output_block(
        log, i,
        period_attr="output_period_bgc_his",
        nrpf_attr="output_nrpf_bgc_his",
        vars_attr="output_variables_bgc"
    )


@register_logparser("output_period_bgc_avg", "output_nrpf_bgc_avg", "output_variables_bgc")
def parse_bgc_avg(log, i):
    """Parse the biogeochemistry averaged output file settings from the logfile."""
    if "bgc :: average file" not in log.lines[i].strip().lower():
        return i
    return _parse_output_block(
        log, i,
        period_attr="output_period_bgc_avg",
        nrpf_attr="output_nrpf_bgc_avg",
        vars_attr="output_variables_bgc"
    )

# ------------------- cstar -------------------

@register_logparser("output_period_cstar", "output_nrpf_cstar", "output_variables_cstar")
def parse_cstar(log, i):
    """Parse the settings associated with the custom cstar_output module from the logfile."""
    if "cstar_output ::" not in log.lines[i]:
        return i
    return _parse_output_block(
        log, i,
        period_attr="output_period_cstar",
        nrpf_attr="output_nrpf_cstar",
        vars_attr="output_variables_cstar"
    )

@register_logparser("cdr_releases")
def parse_cdr_releases(log: "ROMSSimulationLog", i: int) -> int:
    """Parse information on Carbon Dioxide Removal (CDR) perturbations from the logfile."""
    line = log.lines[i]
    if not line.lstrip().startswith("The minimum distance to Release"):
        return i

    if not hasattr(log, "cdr_releases"):
        log.cdr_releases = []

    # we expect the next few lines to contain everything
    # pattern:
    # The minimum distance to Release           N is   <float>
    # This is on rank:         RANK
    # at point          I          J
    # The intended release location was Lon: <lon> Lat: <lat>
    # The release will take place at Lon: <lon> Lat: <lat>

    # release number and distance
    m = re.search(r"Release\s+(\d+)\s+is\s+([0-9.Ee+-]+)", line)
    if not m:
        return i+1
    nr = int(m.group(1))
    dist = float(m.group(2))

    # rank
    rank_line = log.lines[i+1]
    m = re.search(r"rank:\s*(\d+)", rank_line)
    rank = int(m.group(1)) if m else None

    # grid ij
    ij_line = log.lines[i+2]
    m = re.search(r"point\s+(\d+)\s+(\d+)", ij_line)
    ij = (int(m.group(1)), int(m.group(2))) if m else None

    # intended
    intended_line = log.lines[i+3]
    m = re.search(r"Lon:\s*([0-9.Ee+-]+)\s+Lat:\s*([0-9.Ee+-]+)", intended_line)
    intended = (float(m.group(1)), float(m.group(2))) if m else None

    # actual
    actual_line = log.lines[i+4]
    m = re.search(r"Lon:\s*([0-9.Ee+-]+)\s+Lat:\s*([0-9.Ee+-]+)", actual_line)
    actual = (float(m.group(1)), float(m.group(2))) if m else None

    log.cdr_releases.append(
        dict(
            release_number=nr,
            distance=dist,
            grid_location=dict(rank=rank, ij=ij),
            intended_location=dict(lon=intended[0], lat=intended[1]) if intended else None,
            location=dict(lon=actual[0], lat=actual[1]) if actual else None,
        )
    )

    # advance past this 5–6 line block
    return i + 5

def _fix_fort_float(s: str) -> str:
    """Patch Fortran floats missing an exponent marker."""
    # ex: "4.21568551562-03" -> "4.21568551562E-03"
    # only patch if pattern is <digits>.<digits><sign><digits>
    # and there is no "E" already
    if ("E" not in s) and ("e" not in s):
        # last + or - that splits mantissa/exponent occurs after at least 2 digits
        # safe simple: if final +/-### pattern is present, insert 'E'
        m = re.search(r"([0-9])([+-][0-9]+)$", s)
        if m:
            return re.sub(r"([0-9])([+-][0-9]+)$", r"\1E\2", s)
    return s

@register_logparser(
    "time_days", "kinetic_energy", "barotropic_ke",
    "max_adv_cfl", "max_vert_cfl",
    "i_cx", "j_cx", "k_c",
    "timestep_walltime"
)
def parse_step_block(log: "ROMSSimulationLog", i: int) -> int:
    """Parse the timestep-by-timestep output columns (under "STEP") from the logfile."""
    line = log.lines[i].strip()

    # detect the header
    if not (line.startswith("STEP") and "KINETIC_ENRG" in line):
        return i

    # initialise storage
    log.time_days = []
    log.kinetic_energy = []
    log.barotropic_ke = []
    log.max_adv_cfl = []
    log.max_vert_cfl = []
    log.i_cx = []
    log.j_cx = []
    log.k_c = []
    log.timestep_walltime = []

    j = i + 1

    # loop until EOF
    while j < len(log.lines):
        s = log.lines[j].strip()
        parts = s.split()
        parts = [_fix_fort_float(p) for p in parts]
        # candidate line: 10 columns, first is int
        if len(parts) == 10 and parts[0].isdigit():
            try:
                log.time_days.append(float(parts[1]))
                log.kinetic_energy.append(float(parts[2]))
                log.barotropic_ke.append(float(parts[3]))
                log.max_adv_cfl.append(float(parts[4]))
                log.max_vert_cfl.append(float(parts[5]))
                log.i_cx.append(int(parts[6]))
                log.j_cx.append(int(parts[7]))
                log.k_c.append(int(parts[8]))
                log.timestep_walltime.append(float(parts[9]))
            except ValueError:
                pass   # silently ignore malformed rows

        j += 1

    return j  # consumed through EOF


if __name__ == "__main__":
    fpath = Path(sys.argv[1])
    log = ROMSSimulationLog(fpath)
    if log.job_id:
        identifier = log.job_id
    else:
        identifier = datetime.datetime.now().strftime("%Y%m%d%H%M%S")

    if identifier in str(fpath):
        log.to_json(fpath.with_name(fpath.name + ".json"))
    else:
        log.to_json(fpath.parent/f"{fpath.stem}_{identifier}.json")





