import os
import roms_tools as rt
import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import matplotlib.pyplot as plt
import subprocess

from pathlib import Path

import numpy as np
import xarray as xr
import pandas as pd

#--------------------------------------------------------------------------------
# Create upstream input datasets for roms-tools:

def create_roms_tools_inputs(target_dir):
    create_rti_bgc_3d(target_dir)
    create_rti_phys_3d(target_dir)
    create_rti_bgc_surf(target_dir)
    create_rti_phys_surf(target_dir)
    create_rti_tides(target_dir)

def create_rti_bgc_3d(target_dir: Path):
    times = pd.date_range("2010-01-01", periods=2, freq="MS")
    z_t = np.array([500.0, 14500.0], dtype=np.float32)
    lat = np.array([33.5, 35.5], dtype=np.float64)
    lon = np.array([-121.5, -119.5], dtype=np.float64)
    z_t_150m = np.array([500.0, 14500.0], dtype=np.float32)

    coords = {
        "time": times,
        "z_t": z_t,
        "lat": lat,
        "lon": lon,
        "z_t_150m": z_t_150m,
    }

    shape_zt = (len(times), len(z_t), len(lat), len(lon))
    shape_150m = (len(times), len(z_t_150m), len(lat), len(lon))
    dims_zt = ("time", "z_t", "lat", "lon")
    dims_150m = ("time", "z_t_150m", "lat", "lon")

    zt_vars = ["ALK", "ALK_ALT_CO2", "DIC", "DIC_ALT_CO2", "DOC", "DOCr",
               "DON", "DONr", "DOP", "DOPr", "Fe", "Lig", "NH4", "NO3",
               "O2", "PO4", "SiO3"]

    z150m_vars = ["diatC", "diatChl", "diatFe", "diatP", "diatSi",
                  "diazC", "diazChl", "diazFe", "diazP",
                  "spC", "spCaCO3", "spChl", "spFe", "spP", "zooC"]

    data_vars = {}
    for v in zt_vars:
        data_vars[v] = xr.DataArray(np.full(shape_zt, np.nan, dtype=np.float32), dims=dims_zt)
    for v in z150m_vars:
        data_vars[v] = xr.DataArray(np.full(shape_150m, np.nan, dtype=np.float32), dims=dims_150m)

    ds = xr.Dataset(data_vars, coords=coords)

    # z_t variables
    ds["ALK"].values[:, :, 0, :]         = [[[2259.30, 2266.32], [2291.23, 2294.27]],
                                            [[2259.79, 2266.31], [2289.59, 2292.84]]]
    ds["ALK_ALT_CO2"].values[:, :, 0, :] = [[[2259.28, 2266.30], [2291.22, 2294.26]],
                                            [[2259.78, 2266.30], [2289.58, 2292.83]]]
    ds["DIC"].values[:, :, 0, :]      = [[[2026.41, 2028.81], [2176.44, 2181.19]],
                                         [[2028.72, 2032.20], [2171.55, 2177.09]]]
    ds["DIC_ALT_CO2"].values[:, :, 0, :] = [[[1987.60, 1990.27], [2150.13, 2155.50]],
                                            [[1990.68, 1994.49], [2145.26, 2151.49]]]
    ds["DOC"].values[:, :, 0, :]      = [[[38.053, 38.246], [9.553,  9.089]],
                                         [[37.094, 37.096], [10.012, 9.460]]]
    ds["DOCr"].values[:, :, 0, :]     = [[[23.118, 23.060], [22.677, 22.659]],
                                         [[23.097, 23.047], [22.689, 22.670]]]
    ds["DON"].values[:, :, 0, :]      = [[[3.3778, 3.4592], [0.7983, 0.7595]],
                                         [[3.2909, 3.3476], [0.8382, 0.7916]]]
    ds["DONr"].values[:, :, 0, :]     = [[[1.2869, 1.2925], [1.2453, 1.2438]],
                                         [[1.2862, 1.2911], [1.2469, 1.2451]]]
    ds["DOP"].values[:, :, 0, :]      = [[[0.18291, 0.18381], [0.03807, 0.03597]],
                                         [[0.17755, 0.17794], [0.04021, 0.03771]]]
    ds["DOPr"].values[:, :, 0, :]     = [[[0.02542, 0.02545], [0.02425, 0.02420]],
                                         [[0.02541, 0.02543], [0.02429, 0.02424]]]
    ds["Fe"].values[:, :, 0, :]       = [[[0.000489, 0.000554], [0.001137, 0.001427]],
                                         [[0.000501, 0.000586], [0.001139, 0.001441]]]
    ds["Lig"].values[:, :, 0, :]      = [[[0.001204, 0.001208], [0.001829, 0.001839]],
                                         [[0.001216, 0.001226], [0.001815, 0.001828]]]
    ds["NH4"].values[:, :, 0, :]      = [[[0.08197, 0.09139], [0.03600, 0.03866]],
                                         [[0.10656, 0.13054], [0.04011, 0.03945]]]
    ds["NO3"].values[:, :, 0, :]      = [[[0.5848, 0.6530], [15.571, 15.762]],
                                         [[0.8839, 1.0752], [15.225, 15.493]]]
    ds["O2"].values[:, :, 0, :]       = [[[251.52, 248.75], [135.62, 129.69]],
                                         [[254.97, 251.45], [139.18, 132.58]]]
    ds["PO4"].values[:, :, 0, :]      = [[[0.24450, 0.24335], [1.36350, 1.39563]],
                                         [[0.26549, 0.27373], [1.33093, 1.36904]]]
    ds["SiO3"].values[:, :, 0, :]     = [[[6.7365, 6.7148], [24.153, 24.770]],
                                         [[6.9181, 7.0353], [23.479, 24.211]]]

    # z_t_150m variables
    ds["diatC"].values[:, :, 0, :]    = [[[0.5924,  0.6263],  [-0.002546, -0.002034]],
                                         [[0.7186,  0.4739],  [-0.002396, -0.002042]]]
    ds["diatChl"].values[:, :, 0, :]  = [[[0.13961, 0.14527], [-0.000149, -0.000100]],
                                         [[0.16560, 0.10658], [-0.000167, -0.000156]]]
    ds["diatFe"].values[:, :, 0, :]   = [[[1.671e-5, 1.846e-5], [-5.316e-8, -4.015e-8]],
                                         [[1.964e-5, 1.390e-5], [-5.037e-8, -4.108e-8]]]
    ds["diatP"].values[:, :, 0, :]    = [[[0.004305, 0.004546], [-1.843e-5, -1.479e-5]],
                                         [[0.005306, 0.003515], [-1.750e-5, -1.502e-5]]]
    ds["diatSi"].values[:, :, 0, :]   = [[[0.06786, 0.07153], [-0.000282, -0.000227]],
                                         [[0.08383, 0.05588], [-0.000269, -0.000233]]]
    ds["diazC"].values[:, :, 0, :]    = [[[0.026776, 0.028059], [-0.000576, -0.000689]],
                                         [[0.016138, 0.020493], [-0.000610, -0.000723]]]
    ds["diazChl"].values[:, :, 0, :]  = [[[0.003921, 0.004050], [-0.000111, -0.000136]],
                                         [[0.002369, 0.002922], [-0.000116, -0.000142]]]
    ds["diazFe"].values[:, :, 0, :]   = [[[1.575e-6, 1.677e-6], [-3.155e-8, -3.767e-8]],
                                         [[9.508e-7, 1.226e-6], [-3.360e-8, -3.964e-8]]]
    ds["diazP"].values[:, :, 0, :]    = [[[1.929e-4, 2.024e-4], [-4.214e-6, -5.054e-6]],
                                         [[1.170e-4, 1.503e-4], [-4.460e-6, -5.297e-6]]]
    ds["spC"].values[:, :, 0, :]      = [[[0.61605, 0.64067], [-1.265e-3, 9.792e-7]],
                                         [[0.71727, 0.72912], [-6.598e-4, 4.369e-5]]]
    ds["spCaCO3"].values[:, :, 0, :]  = [[[0.025301, 0.027064], [-0.000775, -0.000487]],
                                         [[0.032552, 0.034121], [-0.000601, -0.000421]]]
    ds["spChl"].values[:, :, 0, :]    = [[[0.10992, 0.11303], [-0.000214,  0.000164]],
                                         [[0.12565, 0.12489], [-5.311e-5,  1.558e-4]]]
    ds["spFe"].values[:, :, 0, :]     = [[[1.848e-5, 1.922e-5], [-3.713e-8,  9.120e-10]],
                                         [[2.152e-5, 2.187e-5], [-1.914e-8,  2.062e-9]]]
    ds["spP"].values[:, :, 0, :]      = [[[0.004490, 0.004662], [-1.002e-5,  6.133e-7]],
                                         [[0.005329, 0.005452], [-5.006e-6,  8.930e-7]]]
    ds["zooC"].values[:, :, 0, :]     = [[[1.8614, 1.9732], [0.02186, 0.02520]],
                                         [[1.9863, 2.0369], [0.02571, 0.02775]]]

    ds.to_netcdf(target_dir/"fake_bgc_3d_data.nc")

def create_rti_phys_3d(target_dir: Path):
    depth     = np.array([0.494025, 266.0403], dtype=np.float32)
    latitude  = np.array([34.333332, 34.75], dtype=np.float32)
    longitude = np.array([-120.833336, -120.416664], dtype=np.float32)
    times     = pd.date_range("2010-01-01", periods=2, freq="MS").astype("datetime64[ns]")

    coords = {"time": times, "depth": depth, "latitude": latitude, "longitude": longitude}

    dims_4d = ("time", "depth", "latitude", "longitude")
    dims_3d = ("time", "latitude", "longitude")
    shape_4d = (2, 2, 2, 2)
    shape_3d = (2, 2, 2)

    data_vars = {}
    for v in ["so", "thetao", "uo", "vo"]:
        data_vars[v] = xr.DataArray(np.full(shape_4d, np.nan, dtype=np.float64), dims=dims_4d)
        data_vars["zos"] = xr.DataArray(np.full(shape_3d, np.nan, dtype=np.float64), dims=dims_3d)

    ds = xr.Dataset(data_vars, coords=coords)

    # so
    ds["so"].values[0, 0, :, :] = [[33.24686, 33.21482], [33.20719, np.nan]]
    ds["so"].values[0, 1, :, :] = [[34.33485, 34.17310], [np.nan,   np.nan]]
    ds["so"].values[1, 0, :, :] = [[33.18278, 33.13242], [33.04392, np.nan]]
    ds["so"].values[1, 1, :, :] = [[34.31654, 34.15174], [np.nan,   np.nan]]

    # thetao
    ds["thetao"].values[0, 0, :, :] = [[14.75884, 13.43019], [13.29615, np.nan]]
    ds["thetao"].values[0, 1, :, :] = [[ 8.91174,  8.91467], [np.nan,   np.nan]]
    ds["thetao"].values[1, 0, :, :] = [[14.36259, 13.97732], [14.27250, np.nan]]
    ds["thetao"].values[1, 1, :, :] = [[ 9.06775,  8.95862], [np.nan,   np.nan]]

    # zos (3d - no depth dim)
    ds["zos"].values[0, :, :] = [[0.20844, 0.20753], [0.20569, np.nan]]
    ds["zos"].values[1, :, :] = [[0.16633, 0.17701], [0.18220, np.nan]]

    # uo
    ds["uo"].values[0, 0, :, :] = [[-0.07019,  0.09827], [ 0.04517, np.nan]]
    ds["uo"].values[0, 1, :, :] = [[-0.00671, -0.00244], [np.nan,   np.nan]]
    ds["uo"].values[1, 0, :, :] = [[-0.03784,  0.02380], [-0.04639, np.nan]]
    ds["uo"].values[1, 1, :, :] = [[ 0.00000, -0.00610], [np.nan,   np.nan]]

    # vo
    ds["vo"].values[0, 0, :, :] = [[ 0.14893, -0.09156], [-0.03052, np.nan]]
    ds["vo"].values[0, 1, :, :] = [[ 0.10254,  0.00671], [np.nan,   np.nan]]
    ds["vo"].values[1, 0, :, :] = [[ 0.08057, -0.01648], [-0.01099, np.nan]]
    ds["vo"].values[1, 1, :, :] = [[ 0.08911,  0.00183], [np.nan,   np.nan]]

    ds.to_netcdf(target_dir/"fake_phys_3d_data.nc")

def create_rti_bgc_surf(target_dir: Path):
    times = pd.date_range("2010-01-01", periods=2, freq="MS")
    lat = np.array([33.5, 35.5], dtype=np.float64)
    lon = np.array([-121.5, -119.5], dtype=np.float64)

    coords = {
        "time": times,
        "lat": lat,
        "lon": lon,
        "z_t": np.float32(500.0),
    }

    dims = ("time", "lat", "lon")
    shape = (len(times), len(lat), len(lon))

    data_vars = {}
    for v in ["pCO2SURF", "NHy_FLUX", "NOx_FLUX", "IRON_FLUX", "dust_FLUX_IN"]:
        data_vars[v] = xr.DataArray(np.full(shape, np.nan, dtype=np.float32), dims=dims)

    ds = xr.Dataset(data_vars, coords=coords)
    ds.attrs = {
        "Title": "From the 1x1 degree CESM-MARBL simulation.",
        "original_path": "/glade/campaign/collections/cmip/CMIP6/CESM-HR/FOSI_BGC/L...",
        "regrid_method": "conservative",
    }

    ds["pCO2SURF"].values[:, 0, :]    = [[387.27, 387.54],
                                          [382.06, 383.74]]
    ds["NHy_FLUX"].values[:, 0, :]    = [[1.2647e-12, 2.7269e-12],
                                          [1.1239e-12, 2.4886e-12]]
    ds["NOx_FLUX"].values[:, 0, :]    = [[5.5854e-12, 5.9207e-12],
                                          [4.9581e-12, 5.6438e-12]]
    ds["IRON_FLUX"].values[:, 0, :]   = [[3.3165e-08, 7.3737e-08],
                                          [3.6120e-08, 7.6050e-08]]
    ds["dust_FLUX_IN"].values[:, 0, :] = [[3.1484e-12, 6.6325e-12],
                                           [3.6227e-12, 6.6772e-12]]
    ds.to_netcdf(target_dir/"fake_bgc_surf_data.nc")

def create_rti_phys_surf(target_dir: Path):

    times     = np.array(["2010-01-01", "2010-01-31T23:00:00"], dtype="datetime64[ns]")
    latitude  = np.array([34.75, 34.0], dtype=np.float64)
    longitude = np.array([239.0, 240.25], dtype=np.float64)

    coords = {
        "time":      times,
        "latitude":  latitude,
        "longitude": longitude,
        "expver":    xr.DataArray(["0001", "0001"], dims="time"),
    }

    dims   = ("time", "latitude", "longitude")
    shape  = (2, 2, 2)

    data_vars = {}
    for v in ["tp", "ssr", "strd", "u10", "v10", "d2m", "t2m", "sst"]:
        data_vars[v] = xr.DataArray(np.full(shape, np.nan, dtype=np.float32), dims=dims)

    ds = xr.Dataset(data_vars, coords=coords)
    ds["tp"].values   = np.array([[[0., 0.], [0., 0.]],
                                   [[0., 0.], [0., 0.]]], dtype=np.float32)
    ds["ssr"].values  = np.array([[[537152.,  603008.], [633856.,  656128.]],
                                   [[1533824., 1359232.], [1258048., 1232128.]]], dtype=np.float32)
    ds["strd"].values = np.array([[[1100785.8, 1008657.8], [1088241.8, 1042641.8]],
                                   [[1036525.6,  994477.6], [1085229.5, 1144141.5]]], dtype=np.float32)
    ds["u10"].values  = np.array([[[2.47137,  0.58368], [1.56903, 4.37567]],
                                   [[3.69145,  0.68071], [3.37212, 1.64946]]], dtype=np.float32)
    ds["v10"].values  = np.array([[[-2.93275, -0.56654], [-3.55678, -1.50014]],
                                   [[-4.06624,  0.50896], [-5.14827, -0.22542]]], dtype=np.float32)
    ds["d2m"].values  = np.array([[[281.318, 276.676], [282.010, 284.021]],
                                   [[281.344, 272.317], [282.022, 283.844]]], dtype=np.float32)
    ds["t2m"].values  = np.array([[[285.755, 285.039], [286.417, 286.828]],
                                   [[286.199, 285.072], [286.273, 286.846]]], dtype=np.float32)
    ds["sst"].values  = np.array([[[286.698, np.nan], [287.327, 287.159]],
                                   [[287.048, np.nan], [287.251, 287.491]]], dtype=np.float32)

    ds.to_netcdf(target_dir/"fake_phys_surf_data.nc")


def create_rti_tides(target_dir: Path):
    # Shared grid arrays
    lon_z = np.array([[239.00000675, 239.00000675], [239.83334011, 239.83334011]])
    lat_z = np.array([[34.3333343,   35.00000099],  [34.3333343,   35.00000099]])
    lon_u = np.array([[238.91667342, 238.91667342], [239.75000677, 239.75000677]])
    lat_u = np.array([[34.3333343,   35.00000099],  [34.3333343,   35.00000099]])
    lon_v = np.array([[239.00000675, 239.00000675], [239.83334011, 239.83334011]])
    lat_v = np.array([[34.25000097,  34.91666765],  [34.25000097,  34.91666765]])
    con   = np.array([b'm2  ', b's2  '], dtype="S4")

    tpxo_g = xr.Dataset({
        "mz":    xr.DataArray(np.array([[np.nan, np.nan], [np.nan, 0.0]]), dims=("nx", "ny")),
        "mu":    xr.DataArray(np.array([[np.nan, np.nan], [np.nan, 0.0]]), dims=("nx", "ny")),
        "mv":    xr.DataArray(np.array([[np.nan, np.nan], [np.nan, 0.0]]), dims=("nx", "ny")),
        "lon_z": xr.DataArray(lon_z, dims=("nx", "ny")),
        "lat_z": xr.DataArray(lat_z, dims=("nx", "ny")),
        "lon_u": xr.DataArray(lon_u, dims=("nx", "ny")),
        "lat_u": xr.DataArray(lat_u, dims=("nx", "ny")),
        "lon_v": xr.DataArray(lon_v, dims=("nx", "ny")),
        "lat_v": xr.DataArray(lat_v, dims=("nx", "ny")),
    })

    tpxo_h = xr.Dataset({
        "con":   xr.DataArray(con, dims=("nc",)),
        "lon_z": xr.DataArray(lon_z, dims=("nx", "ny")),
        "lat_z": xr.DataArray(lat_z, dims=("nx", "ny")),
        "hRe":   xr.DataArray(np.array([[[-0.44538313, -0.47057059], [-0.45561337, 0.0]],
                                         [[-0.13853204, -0.14326806], [-0.15103002, 0.0]]]), dims=("nc", "nx", "ny")),
        "hIm":   xr.DataArray(np.array([[[-0.13290478, -0.09993473], [-0.17825195, 0.0]],
                                         [[-0.05957966, -0.04452635], [-0.07763831, 0.0]]]), dims=("nc", "nx", "ny")),
    })

    tpxo_u = xr.Dataset({
        "con":   xr.DataArray(con, dims=("nc",)),
        "lon_u": xr.DataArray(lon_u, dims=("nx", "ny")),
        "lat_u": xr.DataArray(lat_u, dims=("nx", "ny")),
        "lon_v": xr.DataArray(lon_v, dims=("nx", "ny")),
        "lat_v": xr.DataArray(lat_v, dims=("nx", "ny")),
        "URe":   xr.DataArray(np.array([[[15.335916,  5.2130747], [12.584041,  0.0]],
                                         [[10.243252,  2.7856405], [ 7.2664237, 0.0]]]), dims=("nc", "nx", "ny")),
        "UIm":   xr.DataArray(np.array([[[ 3.6844215,  0.53955126], [0.14630099, 0.0]],
                                         [[-0.2055368,  0.10970733], [-1.2532367, 0.0]]]), dims=("nc", "nx", "ny")),
        "VRe":   xr.DataArray(np.array([[[-22.15554,  -10.345664],  [-1.384152,  0.0]],
                                         [[-11.172465,  -5.3363204], [-1.1036928, 0.0]]]), dims=("nc", "nx", "ny")),
        "VIm":   xr.DataArray(np.array([[[-14.589685,  -4.605853],  [-0.27541187, 0.0]],
                                         [[ -2.1566737, -0.793107],  [-0.2983798,  0.0]]]), dims=("nc", "nx", "ny")),
    })

    tpxo_h.to_netcdf(target_dir/"fake_tides_data_h.nc")
    tpxo_u.to_netcdf(target_dir/"fake_tides_data_u.nc")
    tpxo_g.to_netcdf(target_dir/"fake_tides_data_g.nc")

#--------------------------------------------------------------------------------
# Create ROMS inputs using roms-tools

def create_roms_inputs(target_dir: Path):

    create_roms_tools_inputs(target_dir)

    grid = create_roms_grid()
    grid.save(target_dir/"example_input_grid.nc")

    create_roms_physical_boundary_forcing(grid, target_dir)
    create_roms_bgc_boundary_forcing(grid, target_dir)

    create_roms_physical_surface_forcing(grid, target_dir)
    create_roms_bgc_surface_forcing(grid, target_dir)

    create_roms_initial_conditions(grid, target_dir)

    create_roms_tides(grid, target_dir)

    create_roms_cdr_forcing_3d(
        grid,
        irange=range(8,12),
        jrange=range(18,22),
        alk_pert=1e-3,
        target_dir=target_dir
    )

    create_roms_cdr_forcing_dp(
        grid,
        release_lon=-120.64,
        release_lat=34.5349,
        alk_pert=1e-3,
        target_dir=target_dir
    )

    create_roms_cdr_forcing_parm(
        grid,
        release_lon=-120.64,
        release_lat=34.5349,
        release_dep=10,
        alk_pert=1e-3,
        target_dir=target_dir
    )
    river_forcing_ds = create_river_or_pipe_data(
        input_type="river",
        forcing_i=[10,11,12,13,14],
        forcing_j=[20,19,19,18,18],
        grid=grid
    )
    river_forcing_ds.to_netcdf(target_dir/"example_input_river_forcing.nc")
    pipe_forcing_ds = create_river_or_pipe_data(
        input_type="pipe",
        forcing_i=[3,4,5],
        forcing_j=[20,19,18],
        grid=grid
    )
    pipe_forcing_ds.to_netcdf(target_dir/"example_input_pipe_forcing.nc")

    create_flux_forcing_from_bulk(target_dir/"example_input_surface_forcing.nc")

    # Pre-processing
    subprocess.run("for F in example_input_*;do partit 2 2 ${F};done",text=True,shell=True,cwd=target_dir)
    subprocess.run("partit 2 2 cdr_forcing_3d.nc",text=True,shell=True,cwd=target_dir)


def create_roms_grid():
    grid=rt.Grid(
        theta_s=6.0,
        theta_b=6.0,
        hc=25.0,
        N=10,
        nx=39, ny=19,
        center_lon=-120.6475,
        center_lat=34.54,
        rot=-37.4,
        size_x=25,
        size_y=12,
    )
    return grid

def create_roms_physical_boundary_forcing(grid, target_dir):
    rpf=rt.BoundaryForcing(grid=grid,
                   start_time=dt.datetime(2010,1,1),
                   end_time=dt.datetime(2010,1,2),
                   source={"name":"GLORYS","path":target_dir/"fake_phys_3d_data.nc"},
                       apply_2d_horizontal_fill=True
                  )
    rpf.save(target_dir/"example_input_boundary_forcing.nc",group=False)

def create_roms_bgc_boundary_forcing(grid, target_dir):
    rbf=rt.BoundaryForcing(grid=grid,
                           type="bgc",
                           start_time=dt.datetime(2010,1,1),
                           end_time=dt.datetime(2010,1,2),
                           source={"name":"CESM_REGRIDDED",
                                   "path":target_dir/"fake_bgc_3d_data.nc"},
                           apply_2d_horizontal_fill=True
                           )
    # Add DOFE for BEC:
    for bry in ["west","north","south","east"]:
        rbf.ds[f"DOFE_{bry}"] = rbf.ds[f"DOP_{bry}"]
        rbf.ds[f"DOFE_{bry}"][:] = 0.
        rbf.ds[f"DOFE_{bry}"].attrs = {
            "long_name": f"{bry}ern boundary dissolved organic iron",
            "units"    : "mMol Fe m-3",
        }

    rbf.save(target_dir/"example_input_bgc_boundary_forcing.nc",group=False)


def create_roms_physical_surface_forcing(grid, target_dir):
    sf=rt.SurfaceForcing(grid=grid,
                  start_time=dt.datetime(2010,1,1),
                  end_time=dt.datetime(2010,1,2),
                  source={"name":"ERA5","path":target_dir/"fake_phys_surf_data.nc"},
                 )
    sf.save(target_dir/"example_input_surface_forcing.nc",group=False)

def create_roms_bgc_surface_forcing(grid, target_dir):
    rsf=rt.SurfaceForcing(grid=grid,
                          type="bgc",
                          start_time=dt.datetime(2010,1,1),
                          end_time=dt.datetime(2010,1,2),
                          source={"name":"CESM_REGRIDDED",
                                  "path":target_dir/"fake_bgc_surf_data.nc"}
                          )
    rsf.save(target_dir/"example_input_bgc_surface_forcing.nc",group=False)


def create_roms_initial_conditions(grid, target_dir):
    ic = rt.InitialConditions(
        grid=grid,
        ini_time=dt.datetime(2010,1,1),
        source={"name":"GLORYS","path":target_dir/"fake_phys_3d_data.nc"},
        bgc_source={"name":"CESM_REGRIDDED","path":target_dir/"fake_bgc_3d_data.nc"}
    )
    # Add DOFE for BEC
    ic.ds["DOFE"] = ic.ds["DOP"]
    ic.ds["DOFE"][:] = 0.
    ic.ds["DOFE"].attrs = {
        "long_name": "Dissolved organic iron",
        "units": "mMol Fe m-3"
    }


    ic.save(target_dir/"example_input_bgc_initial_conditions.nc")

def create_roms_tides(grid, target_dir):
    tpxo_dict = {
        "grid":  target_dir/"fake_tides_data_g.nc",
        "h":     target_dir/"fake_tides_data_h.nc",
        "u":     target_dir/"fake_tides_data_u.nc"
    }
    tidal_forcing = rt.TidalForcing(
        grid=grid,
        source={"name": "TPXO", "path": tpxo_dict},
        ntides=2,  # Number of constituents to consider <= 15. Default is 10.
        model_reference_date=dt.datetime(2000, 1, 1),  # Model reference date. Default is January 1, 2000.
        use_dask=True
)
    tidal_forcing.save(target_dir/"example_input_tides.nc")


def create_roms_cdr_forcing_3d(grid, irange, jrange, alk_pert, target_dir):

    grid_ds = grid.ds

    cdr_time = xr.DataArray(
            np.array([3653.0, 3655.0]),
            dims=("cdr_time",),
            attrs={"long_name": "Time for CDR release",
                   "units": "days"}
    )
    coords = {
        "cdr_time": cdr_time,
        "s_rho":    grid_ds.s_rho.astype("int32"),
        "xi_rho":   grid_ds.xi_rho.astype("int32"),
        "eta_rho":  grid_ds.eta_rho.astype("int32"),
    }
    shape_3d = (
        len(cdr_time),
        len(grid_ds.s_rho),
        len(grid_ds.eta_rho),
        len(grid_ds.xi_rho)
    )

    alk = xr.DataArray(
        np.zeros(shape_3d, dtype="f8"),
        dims=("cdr_time", "s_rho", "eta_rho", "xi_rho"),
        attrs={
            "long_name": "tracer flux [mmol/s]",
            "units": "mmol/s",
        },
    )
    alk[:, :, irange, jrange] = alk_pert

    dic = xr.DataArray(
        np.zeros(shape_3d, dtype="f8"),
        dims=("cdr_time", "s_rho", "eta_rho", "xi_rho"),
        attrs={
            "long_name": "tracer flux [mmol/s]",
            "units": "mmol/s",
        },
    )
    ds_3d = xr.Dataset(
        data_vars={
            "cdr_trcflx_3d_ALK": alk,
            "cdr_trcflx_3d_DIC": dic,
        },
        coords=coords,
    )
    ds_3d.to_netcdf(target_dir/"cdr_forcing_3d.nc")

def create_roms_cdr_forcing_dp(grid,release_lon,release_lat,alk_pert,target_dir):
    nt = 34
    ncdr = 1

    grid_ds = grid.ds
    cdr_time = xr.DataArray(
            np.array([3653.0, 3655.0]),
            dims=("cdr_time",),
            attrs={"long_name": "Time for CDR release",
                   "units": "days"}
    )

    coords = {
        "s_rho": grid_ds.s_rho,
        "cdr_time": cdr_time,
        "ncdr_prof": ("ncdr_prof", np.arange(1)),
        "nt": ("nt", np.arange(nt)),
        "one": ("one", np.arange(1)),   # scalar dimension used for a string flag
        "two": ("two", np.arange(2)),
    }

    # depth_profiles – a single‑character flag stored as a bytestring (S1)
    depth_profiles = xr.DataArray(
        np.array([b"T"]),
        dims=("one",),
        attrs={"long_name": "depth profiles (T) or Gaussian (F)", "units": "nondim"},
    )
    profile = xr.DataArray(
        np.zeros((len(cdr_time), len(grid_ds.s_rho), 2, ncdr), dtype="f8"),
        dims=("cdr_time", "s_rho", "two", "ncdr_prof"),
        attrs={"long_name": "tracer flux [mmol/s]", "units": "mmol/s"},
    )
    profile[:, :, 0, :] = alk_pert

    # other scalar variables
    cdr_lon = xr.DataArray(
        np.full((ncdr,), release_lon, dtype="f8"),
        dims=("ncdr_prof",),
        attrs={"long_name": "longitude of CDR release [degrees East]", "units": "deg E"},
    )
    cdr_lat = xr.DataArray(
        np.full((ncdr,), release_lat, dtype="f8"),
        dims=("ncdr_prof",),
        attrs={"long_name": "latitude of CDR release [degrees North]", "units": "deg N"},
    )
    cdr_layer_thickness = xr.DataArray(
        np.full((len(cdr_time), len(grid_ds.s_rho), ncdr), 1.0, dtype="f8"),
        dims=("cdr_time", "s_rho", "ncdr_prof"),
        attrs={
            "long_name": "layer thicknesses of CDR release given in a vertical profile [m]",
            "units": "m",
        },
    )

    # Assemble the dataset
    ds_dp = xr.Dataset(
        {
            "depth_profiles": depth_profiles,
            "cdr_trcflx_profile": profile,
            "cdr_lon": cdr_lon,
            "cdr_lat": cdr_lat,
            "cdr_layer_thickness": cdr_layer_thickness,
        },
        coords=coords,
    )
    ds_dp["cdr_time"].attrs["long_name"] = "Time for CDR release"
    ds_dp["cdr_time"].attrs["units"] = "days"
    ds_dp.to_netcdf(target_dir/"cdr_forcing_dp.nc")


def create_roms_cdr_forcing_parm(grid, release_lon, release_lat, release_dep, alk_pert,target_dir):
    ncdr=1
    nt = 34
    grid_ds=grid.ds
    cdr_time = xr.DataArray(
        np.array([3653.0, 3655.0]),
        dims=("cdr_time",),
        attrs={"long_name": "Time for CDR release",
               "units": "days"}
    )

    coords = {
        "cdr_time": cdr_time,
        "ncdr_parm": ("ncdr_parm", np.arange(1)),
        "nt": ("nt", np.arange(nt)),
    }


    trcflx_parm = xr.DataArray(
        np.zeros((len(cdr_time), nt, ncdr), dtype="f8"),
        dims=("cdr_time", "nt", "ncdr_parm"),
        attrs={"long_name": "tracer flux [mmol/s]", "units": "mmol/s"},
    )
    trcflx_parm[:, 11, :] = alk_pert

    # other scalar variables
    cdr_lon = xr.DataArray(
        np.full((ncdr,), release_lon, dtype="f8"),
        dims=("ncdr_prof",),
        attrs={"long_name": "longitude of CDR release [degrees East]", "units": "deg E"},
    )
    cdr_lat = xr.DataArray(
        np.full((ncdr,), release_lat, dtype="f8"),
        dims=("ncdr_prof",),
        attrs={"long_name": "latitude of CDR release [degrees North]", "units": "deg N"},
    )
    cdr_dep = xr.DataArray(
        np.full((ncdr,), release_dep, dtype="f8"),
        dims=("ncdr_parm",),
        attrs={"long_name": "depth of CDR release [m]", "units": "m"},
    )

    cdr_hsc = xr.DataArray(
        np.full((ncdr,), 360.0, dtype="f8"),
        dims=("ncdr_parm",),
        attrs={"long_name": "horizontal scale of CDR release", "units": "m"},
    )

    cdr_vsc = xr.DataArray(
        np.full((ncdr,), 5.0, dtype="f8"),
        dims=("ncdr_parm",),
        attrs={"long_name": "vertical scale of CDR release", "units": "m"},
    )

    # Assemble
    ds_parm = xr.Dataset(
        {
            "cdr_trcflx": trcflx_parm,
            "cdr_lon"   : cdr_lon,
            "cdr_lat"   : cdr_lat,
            "cdr_dep"   : cdr_dep,
            "cdr_hsc"   : cdr_hsc,
            "cdr_vsc"   : cdr_vsc,
        },
        coords = coords
    )

    ds_parm.to_netcdf(target_dir/"cdr_forcing_parm.nc")

def create_river_or_pipe_data(
    input_type,
    forcing_i,
    forcing_j,
    grid
):
    grid_ds = grid.ds
    n_inputs = 1
    ntracers = 2
    forcing_time=np.arange(0.5, 360.5)

    # Time handling
    ref_time = pd.Timestamp("2000-01-01 00:00:00")
    abs_time = ref_time + pd.to_timedelta(forcing_time, unit="D")

    # Common coordinates
    forcing_time_coords = {
        f"{input_type}_time": forcing_time,
        "abs_time": (f"{input_type}_time", abs_time),
        f"n{input_type}": np.arange(1, n_inputs + 1),
        f"{input_type}_name": (
            f"n{input_type}",
            [f"Example{input_type.capitalize()}"],
        ),
    }
    tracer_coords = {
        "ntracers": np.arange(1, ntracers + 1),
        "tracer_name": ("ntracers", ["temp", "salt"]),
        "tracer_unit": ("ntracers", ["degrees Celsius", "PSU"]),
        "tracer_long_name": ("ntracers", ["potential temperature", "salinity"]),
    }

    # Forcing volume
    constant_value = 2000.  # m3/s
    forcing_volume_data = np.full((len(forcing_time), n_inputs), constant_value, dtype=np.float32)
    forcing_volume = xr.DataArray(
        forcing_volume_data,
        dims=(f"{input_type}_time", f"n{input_type}"),
        coords=forcing_time_coords,
        attrs={
            "long_name": f"{input_type.capitalize()} volume flux",
            "units": "m^3/s",
        },
    )
    forcing_volume[f"{input_type}_time"].attrs.update(
        {
            "long_name": "relative time: days since 2000-01-01 00:00:00",
            "units": "yearday",
            "cycle_length": 360
        }
    )

    # Forcing tracer
    forcing_tracer_data = np.zeros(
        (len(forcing_time), ntracers, n_inputs), dtype=np.float32
    )
    forcing_tracer_data[:, 0, :] = 17.0  # Temperature
    forcing_tracer_data[:, 1, :] = 1.0  # Salinity
    forcing_tracer = xr.DataArray(
        forcing_tracer_data,
        dims=(f"{input_type}_time", "ntracers", f"n{input_type}"),
        coords={**forcing_time_coords,**tracer_coords},
        attrs={
            "long_name": f"{input_type.capitalize()} tracer data",
        },
    )
    forcing_tracer[f"{input_type}_time"].attrs.update(forcing_volume[f"{input_type}_time"].attrs)

    # Forcing index
    forcing_index_data = np.zeros((len(grid_ds.eta_rho), len(grid_ds.xi_rho)))
    forcing_index_data[forcing_i, forcing_j] = 1.
    forcing_index = xr.DataArray(
        forcing_index_data,
        dims=("eta_rho", "xi_rho"),
        attrs={
            "long_name": f"{input_type.capitalize()} ID",
            "units": "none",
        },
    )

    # Forcing fraction
    forcing_fraction_data = np.zeros((len(grid_ds.eta_rho), len(grid_ds.xi_rho)))
    forcing_fraction_data[forcing_i, forcing_j] = 1./len(forcing_i)
    forcing_fraction = xr.DataArray(
        forcing_fraction_data,
        dims=("eta_rho", "xi_rho"),
        attrs={
            "long_name": f"{input_type.capitalize()} volume fraction",
            "units": "none",
        },
    )

    # Assemble Dataset
    forcing_ds = xr.Dataset(
        data_vars={
            f"{input_type}_volume": forcing_volume,
            f"{input_type}_tracer": forcing_tracer,
            f"{input_type}_index": forcing_index,
            f"{input_type}_fraction": forcing_fraction,
        },
        attrs={
        },
    )

    return forcing_ds


def create_flux_forcing_from_bulk(bulk_file):


    bulk_ds = xr.open_dataset(bulk_file, decode_times=False)
    # ---------------------------
    # helpers: rho → u/v
    # ---------------------------
    def rho_to_u(A):
        return (
            0.5 * (A[..., :-1] + A[..., 1:])
        ).rename({"xi_rho": "xi_u"})

    def rho_to_v(A):
        return (
            0.5 * (A[:, :-1, :] + A[:, 1:, :])
        ).rename({"eta_rho": "eta_v"})

    # ---------------------------
    # stresses (new arrays)
    # ---------------------------
    rho_air = 1.22
    Cd = 1.3e-3

    U = bulk_ds["uwnd"]
    V = bulk_ds["vwnd"]
    W = np.hypot(U, V)

    tau_x = rho_air * Cd * W * U
    tau_y = rho_air * Cd * W * V

    sustr = rho_to_u(tau_x).astype("float32")
    svstr = rho_to_v(tau_y).astype("float32")

    # ---------------------------
    # fluxes (new arrays)
    # ---------------------------
    shflux = (bulk_ds["swrad"] + bulk_ds["lwrad"]).astype("float32")
    swrad  = bulk_ds["swrad"].astype("float32")

    # freshwater flux: E − P (placeholder)
    swflux = (-bulk_ds["rain"]).astype("float32")

    # ---------------------------
    # SST / SSS / dQdSST (synthetic)
    # ---------------------------
    template = bulk_ds["Tair"]

    SST     = xr.full_like(template, 20.0).astype("float32")
    SSS     = xr.full_like(template, 35.0).astype("float32")
    dQdSST  = xr.full_like(template, -30.0).astype("float32")

    # ---------------------------
    # time
    # ---------------------------
    time = bulk_ds["time"].astype("float32")

    # ---------------------------
    # New dataset
    # ---------------------------
    forcing = xr.Dataset(
        data_vars=dict(
            sustr   = (("sms_time","eta_rho","xi_u"),   sustr.data),
            svstr   = (("sms_time","eta_v","xi_rho"),   svstr.data),
            shflux  = (("shf_time","eta_rho","xi_rho"), shflux.data),
            swflux  = (("swf_time","eta_rho","xi_rho"), swflux.data),
            SST     = (("sst_time","eta_rho","xi_rho"), SST.data),
            SSS     = (("sss_time","eta_rho","xi_rho"), SSS.data),
            dQdSST  = (("sst_time","eta_rho","xi_rho"), dQdSST.data),
            swrad   = (("srf_time","eta_rho","xi_rho"), swrad.data),
        ),
        coords=dict(
            sms_time = ("sms_time", time.data),
            shf_time = ("shf_time", time.data),
            swf_time = ("swf_time", time.data),
            sst_time = ("sst_time", time.data),
            sss_time = ("sss_time", time.data),
            srf_time = ("srf_time", time.data),
        ),
    )

    forcing.to_netcdf(bulk_file.parent/"example_input_surface_flux_forcing.nc")

create_roms_inputs(Path("input_data"))
