import netCDF4 as netcdf
import numpy as np

def create_netcdf(filename,max_count):
  fh2 = netcdf.Dataset(filename,'w',format='NETCDF4')
  fh2.createDimension('Nobs',None)
  x1 = fh2.createVariable('n', np.int32, ('Nobs'))
  x1.long_name     = 'Profile'
  x1.standard_name = 'no'
  
  fh2.createDimension('pressure',2000)
  x1 = fh2.createVariable('pressure', 'f4', ('pressure'))
  x1.long_name     = 'Pressure'
  x1.standard_name = 'depth'
  x1.units         = 'dbar'
  x1[:]            = np.arange(2000)

  fh2.createDimension('mixed_layers',max_count)
  x1 = fh2.createVariable('mixed_layers', 'f4', ('mixed_layers'),fill_value=0)
  x1.long_name     = 'Mixed-Layer-count'
  x1.standard_name = 'MLD'
  x1[:]            = np.arange(0,np.int(max_count))

  fh2.createDimension('gradient_layers',max_count-1)
  x1 = fh2.createVariable('gradient_layers', 'f4', ('gradient_layers'),fill_value=0)
  x1.long_name     = 'Gradient-Layer-count'
  x1.standard_name = 'MLD'
  x1[:]            = np.arange(0,np.int(max_count)-1)
  
  x2 = fh2.createVariable('prof',np.int32, ('Nobs'),fill_value=0)
  x2.long_name     = 'Profile number of float'
  x2.standard_name = 'prof'
  
  x2 = fh2.createVariable('FloatID',np.int32, ('Nobs'),fill_value=0)
  x2.long_name     = 'Float ID'
  x2.standard_name = 'FloatID'
  
  x2 = fh2.createVariable('juld','f4',('Nobs'),fill_value=0)
  x2.long_name     = 'Julian date of profile'
  x2.standard_name = 'juld'
  x2.units         = 'days after 1950-01-01'
  
  x2 = fh2.createVariable('lon','f4',('Nobs'),fill_value=0)
  x2.long_name     = 'Longitude of float'
  x2.standard_name = 'lon'
  x2.units         = 'degrees'
  
  x2 = fh2.createVariable('lat','f4',('Nobs'),fill_value=0)
  x2.long_name     = 'Latitude of float'
  x2.standard_name = 'lat'
  x2.units         = 'degrees'

  # other variables
  x2 = fh2.createVariable('ct','f4',('Nobs','pressure'),fill_value=0)
  x2.long_name     = 'Conservative Temperature'
  x2.standard_name = 'ct'
  x2.units         = 'degrees Celsius'

  x2 = fh2.createVariable('sa','f4',('Nobs','pressure'),fill_value=0)
  x2.long_name     = 'Absolute Salinity'
  x2.standard_name = 'sa'
  x2.units         = 'g/kg'

  # mixed layers
  x2 = fh2.createVariable('ml_T','f4',('Nobs','mixed_layers'),fill_value=0)
  x2.long_name     = 'average temperature of the mixed layer'
  x2.standard_name = 'ml_T'
  x2.units         = 'degrees Celsius'
 
  x2 = fh2.createVariable('ml_S','f4',('Nobs','mixed_layers'),fill_value=0)
  x2.long_name     = 'average salinity of the mixed layer'
  x2.standard_name = 'ml_S'
  x2.units         = 'g kg-1'
 
  x2 = fh2.createVariable('ml_r','f4',('Nobs','mixed_layers'),fill_value=0)
  x2.long_name     = 'average density of the mixed layer (sigma1)'
  x2.standard_name = 'ml_r'
  x2.units         = 'kg m-3'

  x2 = fh2.createVariable('ml_p','f4',('Nobs','mixed_layers'),fill_value=0)
  x2.long_name     = 'average pressure of the mixed layer'
  x2.standard_name = 'ml_p'
  x2.units         = 'dbar'

  x2 = fh2.createVariable('ml_dT','f4',('Nobs','mixed_layers'),fill_value=0)
  x2.long_name     = 'temperature difference within the mixed layer'
  x2.standard_name = 'ml_dT'
  x2.units         = 'degrees Celsius'

  x2 = fh2.createVariable('ml_dS','f4',('Nobs','mixed_layers'),fill_value=0)
  x2.long_name     = 'salinity difference within the mixed layer'
  x2.standard_name = 'ml_dS'
  x2.units         = 'g kg-1'

  x2 = fh2.createVariable('ml_dr','f4',('Nobs','mixed_layers'),fill_value=0)
  x2.long_name     = 'density difference within the mixed layer'
  x2.standard_name = 'ml_dr'
  x2.units         = 'kg m-3'

  x2 = fh2.createVariable('ml_h','f4',('Nobs','mixed_layers'),fill_value=0)
  x2.long_name     = 'height of the mixed layer'
  x2.standard_name = 'ml_h'
  x2.units         = 'dbar'
 
  x2 = fh2.createVariable('ml_Tu','f4',('Nobs','mixed_layers'),fill_value=0)
  x2.long_name     = 'average Turner angle of the mixed layer'
  x2.standard_name = 'ml_Tu'
  x2.units         = 'degrees'

  x2 = fh2.createVariable('ml_R','f4',('Nobs','mixed_layers'),fill_value=0)
  x2.long_name     = 'average density ratio of the mixed layer'
  x2.standard_name = 'ml_R'
  x2.units         = ' '

  x2 = fh2.createVariable('ml_count',np.int32,('Nobs','mixed_layers'))
  x2.long_name     = 'number of identified gradient layer (low = top, high = deep)'
  x2.standard_name = 'ml_count'
  x2.units         = ' '


  #gradient layer variables
  x2 = fh2.createVariable('gl_dT','f4',('Nobs','gradient_layers'),fill_value=0)
  x2.long_name     = 'temperature difference in gradient layer'
  x2.standard_name = 'gl_dT'
  x2.units         = 'degrees Celsius'

  x2 = fh2.createVariable('gl_dS','f4',('Nobs','gradient_layers'),fill_value=0)
  x2.long_name     = 'salinity difference in gradient layer'
  x2.standard_name = 'gl_dS'
  x2.units         = 'g kg-1'

  x2 = fh2.createVariable('gl_dr','f4',('Nobs','gradient_layers'),fill_value=0)
  x2.long_name     = 'density difference in gradient layer (sigma1)'
  x2.standard_name = 'gl_dr'
  x2.units         = 'kg m-3 dbar-1'

  x2 = fh2.createVariable('gl_h','f4',('Nobs','gradient_layers'),fill_value=0)
  x2.long_name     = 'gradient layer height'
  x2.standard_name = 'gl_h'
  x2.units         = 'dbar'

  x2 = fh2.createVariable('gl_Tu','f4',('Nobs','gradient_layers'),fill_value=0)
  x2.long_name     = 'average Turner angle of the gradient layer'
  x2.standard_name = 'gl_Tu'
  x2.units         = 'degrees'

  x2 = fh2.createVariable('gl_R','f4',('Nobs','gradient_layers'),fill_value=0)
  x2.long_name     = 'average density ratio of the gradient layer'
  x2.standard_name = 'gl_R'
  x2.units         = ' '

  x2 = fh2.createVariable('gl_count',np.int32,('Nobs','gradient_layers'),fill_value=0)
  x2.long_name     = 'number of identified gradient layer (low = top, high = deep)'
  x2.standard_name = 'gl_count'
  x2.units         = ' '

  x2 = fh2.createVariable('gl_dTdz','f4',('Nobs','gradient_layers'),fill_value=0)
  x2.long_name     = 'temperature gradient in gradient layer'
  x2.standard_name = 'gl_dTdz'
  x2.units         = 'degrees Celsius dbar-1'

  x2 = fh2.createVariable('gl_dSdz','f4',('Nobs','gradient_layers'),fill_value=0)
  x2.long_name     = 'temperature difference in gradient layer'
  x2.standard_name = 'gl_dSdz'
  x2.units         = 'g kg-1 dbar-1'


  # mask variables
  x2 = fh2.createVariable('mask_gl_num',np.int32,('Nobs','pressure'),fill_value=0)
  x2.long_name     = 'mask with all gradient layers'
  x2.standard_name = 'mask_gl_num'
  x2.units         = ' '

  x2 = fh2.createVariable('mask_ml_num',np.int32,('Nobs','pressure'),fill_value=0)
  x2.long_name     = 'mask with all mixed layers'
  x2.standard_name = 'mask_ml_num'
  x2.units         = ' '

  x2 = fh2.createVariable('mask_regime',np.int32,('Nobs','gradient_layers'),fill_value=0)
  x2.long_name     = 'regime of gradient layer (DC = 2; SF = -2)'
  x2.standard_name = 'mask_regime'
  x2.units         = ' '

  x2 = fh2.createVariable('mask_qc',np.int32,('Nobs','gradient_layers'),fill_value=0)
  x2.long_name     = 'quality check of gradient layers'
  x2.standard_name = 'mask_qc'
  x2.units         = ' '

  x2 = fh2.createVariable('mask_seq_sf_num',np.int32,('Nobs','gradient_layers'),fill_value=0)
  x2.long_name     = 'number of the sequence in salt-finger regime'
  x2.standard_name = 'mask_seq_sf_num'
  x2.units         = ' '

  x2 = fh2.createVariable('mask_seq_sf_layer',np.int32,('Nobs','gradient_layers'),fill_value=0)
  x2.long_name     = 'layer number in sequence in salt-finger regime'
  x2.standard_name = 'mask_seq_sf_layer'
  x2.units         = ' '

  x2 = fh2.createVariable('mask_seq_sf_count',np.int32,('Nobs','gradient_layers'),fill_value=0)
  x2.long_name     = 'total number of layers per sequence in salt-finger regime'
  x2.standard_name = 'mask_seq_sf_count'
  x2.units         = ' '

  x2 = fh2.createVariable('mask_seq_dc_num',np.int32,('Nobs','gradient_layers'),fill_value=0)
  x2.long_name     = 'number of the sequence in diffusive-convection regime'
  x2.standard_name = 'mask_dc_num'
  x2.units         = ' '

  x2 = fh2.createVariable('mask_seq_dc_layer',np.int32,('Nobs','gradient_layers'),fill_value=0)
  x2.long_name     = 'layer number in sequence in  diffusive-convection regime'
  x2.standard_name = 'mask_dc_layer'
  x2.units         = ' '

  x2 = fh2.createVariable('mask_seq_dc_count',np.int32,('Nobs','gradient_layers'),fill_value=0)
  x2.long_name     = 'total number of layers per sequence in diffusive-convection regime'
  x2.standard_name = 'mask_dc_count'
  x2.units         = ' '

  x2 = fh2.createVariable('mask_gl_sf_layer',np.int32,('Nobs','gradient_layers'),fill_value=0)
  x2.long_name     = 'mask with sequences of gradient layers in salt-finger regime'
  x2.standard_name = 'mask_gl_sf_layer'
  x2.units         = ' '

  x2 = fh2.createVariable('mask_ml_sf_layer',np.int32,('Nobs','mixed_layers'),fill_value=0)
  x2.long_name     = 'mask with sequences of mixed layers in salt-finger regime'
  x2.standard_name = 'mask_ml_sf_layer'
  x2.units         = ' '

  x2 = fh2.createVariable('mask_gl_sf',np.int32,('Nobs','pressure'),fill_value=0)
  x2.long_name     = 'big mask with sequences of gradient layers in salt-finger regime'
  x2.standard_name = 'mask_gl_sf'
  x2.units         = ' '

  x2 = fh2.createVariable('mask_ml_sf',np.int32,('Nobs','pressure'),fill_value=0)
  x2.long_name     = 'big mask with sequences of mixed layers in salt-finger regime'
  x2.standard_name = 'mask_ml_sf'
  x2.units         = ' '

  x2 = fh2.createVariable('mask_gl_dc_layer',np.int32,('Nobs','gradient_layers'),fill_value=0)
  x2.long_name     = 'mask with sequences of gradient layers in diffusive-convection regime'
  x2.standard_name = 'mask_gl_dc_layer'
  x2.units         = ' '


  x2 = fh2.createVariable('mask_ml_dc_layer',np.int32,('Nobs','mixed_layers'),fill_value=0)
  x2.long_name     = 'mask with sequences of mixed layers in diffusive-convection regime'
  x2.standard_name = 'mask_ml_dc_layer'
  x2.units         = ' '


  x2 = fh2.createVariable('mask_gl_dc',np.int32,('Nobs','pressure'),fill_value=0)
  x2.long_name     = 'big mask with sequences of gradient layers in diffusive-convection regime'
  x2.standard_name = 'mask_gl_dc'
  x2.units         = ' '

  x2 = fh2.createVariable('mask_ml_dc',np.int32,('Nobs','pressure'),fill_value=0)
  x2.long_name     = 'big mask with sequences of mixed layers in diffusive-convection regime'
  x2.standard_name = 'mask_ml_dc'
  x2.units         = ' '

  fh2.close()
