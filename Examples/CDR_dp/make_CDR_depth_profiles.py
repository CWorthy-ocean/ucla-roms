import netCDF4

s_rho = 50
times = 2
ncdr = 1
nt = 34

f = netCDF4.Dataset('cdr_forcing_dp.nc', 'w', format='NETCDF4')

tmp = f.createDimension('s_rho', s_rho)
tmp = f.createDimension('cdr_time', times)
tmp = f.createDimension('ncdr_prof', ncdr)
tmp = f.createDimension('nt', nt)
tmp = f.createDimension('one', 1)
tmp = f.createDimension('two', 2)

tmp_zi = f.createVariable('depth_profiles', 'S1', ('one'))
setattr(tmp_zi, 'long_name', "depth profiles (T) or Gaussian (F)")
setattr(tmp_zi, 'units', "nondim")
f['depth_profiles'][:] = 'T'

tmp_zi = f.createVariable('cdr_trcflx_profile', 'f8', ('cdr_time','s_rho','two','ncdr_prof'))
setattr(tmp_zi, 'long_name', "tracer flux [mmol/s]")
setattr(tmp_zi, 'units', "mmol/s")
f['cdr_trcflx_profile'][:,:,:,:] = 0
f['cdr_trcflx_profile'][:,:,0,:] = 1e-3

tmp_zi = f.createVariable('cdr_time', 'f8', ('cdr_time'))
setattr(tmp_zi, 'long_name', "Time for CDR release")
setattr(tmp_zi, 'units', "days")
f['cdr_time'][0] = 3653.0
f['cdr_time'][1] = 3655.0

tmp_zi = f.createVariable('cdr_lon', 'f8', ('ncdr_prof'))
setattr(tmp_zi, 'long_name', "longitude of CDR release [degrees East]")
setattr(tmp_zi, 'units', "deg E")
f['cdr_lon'][:] = -120.64

tmp_zi = f.createVariable('cdr_lat', 'f8', ('ncdr_prof'))
setattr(tmp_zi, 'long_name', "latitude of CDR release [degrees North]")
setattr(tmp_zi, 'units', "deg N")
f['cdr_lat'][:] = 34.5349

tmp_zi = f.createVariable('cdr_layer_thickness', 'f8', ('cdr_time','s_rho','ncdr_prof'))
setattr(tmp_zi, 'long_name', "layer thicknesses of CDR release given in a vertical profile [m]")
setattr(tmp_zi, 'units', "m")
f['cdr_layer_thickness'][:,:,:] = 1

f.close()
