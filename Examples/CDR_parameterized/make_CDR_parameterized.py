import netCDF4

s_rho = 50
ncdr = 1
nt = 34
times = 2

f = netCDF4.Dataset('cdr_forcing_parm.nc', 'w', format='NETCDF4')

tmp = f.createDimension('cdr_time', times)
tmp = f.createDimension('ncdr_parm', ncdr)
tmp = f.createDimension('nt', nt)

tmp_zi = f.createVariable('cdr_trcflx', 'f8', ('cdr_time','nt','ncdr_parm'))
setattr(tmp_zi, 'long_name', "tracer flux [mmol/s]")
setattr(tmp_zi, 'units', "mmol/s")
f['cdr_trcflx'][:,:,:] = 0
f['cdr_trcflx'][:,11,:] = 1e-3

tmp_zi = f.createVariable('cdr_time', 'f8', ('cdr_time'))
setattr(tmp_zi, 'long_name', "Time for CDR release")
setattr(tmp_zi, 'units', "days")
f['cdr_time'][0] = 3653.0
f['cdr_time'][1] = 3655.0

tmp_zi = f.createVariable('cdr_lon', 'f8', ('ncdr_parm'))
setattr(tmp_zi, 'long_name', "longitude of CDR release [degrees East]")
setattr(tmp_zi, 'units', "deg E")
f['cdr_lon'][:] = -120.64

tmp_zi = f.createVariable('cdr_lat', 'f8', ('ncdr_parm'))
setattr(tmp_zi, 'long_name', "latitude of CDR release [degrees North]")
setattr(tmp_zi, 'units', "deg E")
f['cdr_lat'][:] = 34.5349

tmp_zi = f.createVariable('cdr_dep', 'f8', ('ncdr_parm'))
setattr(tmp_zi, 'long_name', "depth of CDR release [m]")
setattr(tmp_zi, 'units', "m")
f['cdr_dep'][:] = 10.0

tmp_zi = f.createVariable('cdr_hsc', 'f8', ('ncdr_parm'))
setattr(tmp_zi, 'long_name', "horizontal scale of CDR release")
setattr(tmp_zi, 'units', "m")
f['cdr_hsc'][:] = 360.0

tmp_zi = f.createVariable('cdr_vsc', 'f8', ('ncdr_parm'))
setattr(tmp_zi, 'long_name', "vertical scale of CDR release")
setattr(tmp_zi, 'units', "m")
f['cdr_vsc'][:] = 5.0

f.close()
