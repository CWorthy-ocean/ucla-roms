import netCDF4
s_rho = 50
xi_rho = 201
eta_rho = 101
times = 2

f = netCDF4.Dataset('cdr_forcing_3d.nc', 'w', format='NETCDF4')

tmp = f.createDimension('xi_rho', xi_rho)
tmp = f.createDimension('eta_rho', eta_rho)
tmp = f.createDimension('s_rho', s_rho)
tmp = f.createDimension('cdr_time', times)

tmp_zi = f.createVariable('cdr_trcflx_3d_ALK', 'f8', ('cdr_time','s_rho','eta_rho','xi_rho'))
setattr(tmp_zi, 'long_name', "tracer flux [mmol/s]")
setattr(tmp_zi, 'units', "mmol/s")
f['cdr_trcflx_3d_ALK'][:,:,:,:] = 0
f['cdr_trcflx_3d_ALK'][:,:,40:60,90:110] = 1e-3

tmp_zi = f.createVariable('cdr_trcflx_3d_DIC', 'f8', ('cdr_time','s_rho','eta_rho','xi_rho'))
setattr(tmp_zi, 'long_name', "tracer flux [mmol/s]")
setattr(tmp_zi, 'units', "mmol/s")
f['cdr_trcflx_3d_DIC'][:,:,:,:] = 0

tmp_zi = f.createVariable('cdr_time', 'f8', ('cdr_time'))
setattr(tmp_zi, 'long_name', "Time for CDR release")
setattr(tmp_zi, 'units', "days")
f['cdr_time'][0] = 3653.0
f['cdr_time'][1] = 3655.0

f.close()
