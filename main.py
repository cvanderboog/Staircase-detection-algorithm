import netCDF4 as netcdf
import numpy as np
import scipy.interpolate as interpolate
import scipy.stats as stats
import scipy.ndimage as ndimage
import os
import shutil
import zipfile

#import functions
from data_preparation import *
from create_netcdf import *
from staircase_detection_algorithm import *

"""
Script to download profiles from ftp server, to interpolate them to a 1 dbar resolution and to detect 
thermohaline staircases in the interpolated profiles. A detailed description can be found in:

van der Boog, C.G. et al. (20xx), Global dataset of thermohaline staircases obtained from Argo
floats and Ice Tethered Profilers. Submitted to Earth System Science Data

made by: Carine van der Boog (2020)
"""

for flt_type in range(1):
  copyright()
  if flt_type == 0:
    print('Argo float profiles')
    centers = ['aoml','bodc','coriolis','csio','csiro','incois','jma','kma','kordi','meds','nmdis']
    directory,floats,list1 = data_list_argo(centers,'Floatlist_argo.xlsx')
  else:
    print('Ice tethered profiles')
    list1 = data_list_itp()
    floats=list1 
 
  for i in range(len(list1)):
    if flt_type == 0:
      if i in np.arange(0,100000,250):
        #create netcdf file to save variables (250 floats per netcdf)
        ncfile = 'argo_'+str(i).zfill(8)+'_'+str(np.min([i+249,len(list1)])).zfill(8)+'.nc'
        create_netcdf(ncfile,200)
    else:
      if i in np.arange(0,100000,10):
        #create netcdf file to save variables (10 floats per netcdf)
        ncfile = 'itp_'+str(i).zfill(8)+'_'+str(np.min([i+9,len(list1)])).zfill(8)+'.nc'
        create_netcdf(ncfile,200)
  
    #download data from ftp server and interpolate to vertical resolution of 1 dbar
    if flt_type == 0:
      if len(np.intersect1d(floats,list1[i]))>0:
        index = np.where(floats==list1[i])[0][0]
        docname1 = 'ftp://ftp.ifremer.fr/ifremer/argo/dac/'+centers[directory[index]]+'/'+str(list1[i])+'/'+str(list1[i])+'_prof.nc'
        filename1 = str(list1[i])+'_prof.nc'
        wget.download(docname1)
        p,lat,lon,ct,sa,juld,z,rho,alpha,beta,sigma = load_data(filename1,True)
        prof_no = int(float(list1[i]))*np.ones(len(lat))
        os.remove(str(list1[i])+'_prof.nc')
    else:
      index = np.where(floats==list1[i])[0][0]
      docname1 = 'ftp://ftp.whoi.edu/whoinet/itpdata/'+list1[i]
      filename1 = list1[i]
      wget.download(docname1)
          
      os.makedirs('tmp/')
      zip_ref = zipfile.ZipFile(floats[index],'r')
      zip_ref.extractall('tmp/')  
      os.remove(floats[index])
      
      profiles = os.listdir('tmp/')
      prof_no,p,lat,lon,ct,sa,juld,z,rho,alpha,beta,sigma = load_data_itp('tmp/',profiles,True)
      shutil.rmtree('tmp/')  
    
    if len(lat)>1:
      #remove profiles with only nans
      n     = np.arange(len(lat),dtype=np.int32)
      n     = n[~np.all(np.isnan(ct), axis=1)]
      p     = p[~np.all(np.isnan(ct), axis=1),:]
      lon   = lon[~np.all(np.isnan(ct), axis=1)]
      lat   = lat[~np.all(np.isnan(ct), axis=1)]
      sa    = sa[~np.all(np.isnan(ct), axis=1),:]
      juld  = juld[~np.all(np.isnan(ct), axis=1)]
      z     = z[~np.all(np.isnan(ct), axis=1),:]
      rho   = rho[~np.all(np.isnan(ct), axis=1),:]
      sigma = sigma[~np.all(np.isnan(ct), axis=1),:]
      alpha = alpha[~np.all(np.isnan(ct), axis=1),:]
      beta  = beta[~np.all(np.isnan(ct), axis=1),:]
      prof_no = prof_no[~np.all(np.isnan(ct), axis=1)]
      ct    = ct[~np.all(np.isnan(ct), axis=1),:]

    if len(p)>0:
      if len(lat)==1 and p.max()>0:
        ct = np.ma.squeeze(ct)[np.newaxis,:]
        sa = np.ma.squeeze(sa)[np.newaxis,:]   
        p  = np.ma.squeeze(p)[np.newaxis,:]
      
      if p.max()>0:
        #detection algorithm
        c1 = 0.0005
        c2 = 0.005
        c3 = 200
        c4 = 30
        gl, ml, masks = get_mixed_layers(np.ma.copy(p),np.ma.copy(ct),np.ma.copy(sa),c1,c2,c3,c4) 
  
        fh2 = netcdf.Dataset(ncfile,'r+')
        t0 = len(fh2.variables['n'][:])
        t1 = len(fh2.variables['n'][:])+len(lat)
       
        #general
        fh2.variables['n'][t0:t1]                   = np.arange(len(lat),dtype=np.int32)
        fh2.variables['lat'][t0:t1]                 = lat
        fh2.variables['lon'][t0:t1]                 = lon
        fh2.variables['prof'][t0:t1]                = np.arange(len(lat))
        fh2.variables['juld'][t0:t1]                = juld
        fh2.variables['ct'][t0:t1,:]                = ct
        fh2.variables['sa'][t0:t1,:]                = sa
        fh2.variables['FloatID'][t0:t1]             = prof_no
        
        #masks
        fh2.variables['mask_gl_num'][t0:t1,:]       = masks.gl_num

        fh2.variables['mask_ml_num'][t0:t1,:]       = masks.ml_num
        fh2.variables['mask_regime'][t0:t1,:]       = masks.regime
        fh2.variables['mask_qc'][t0:t1,:]           = masks.qc
        fh2.variables['mask_seq_sf_num'][t0:t1,:]   = masks.seq_sf_num
        fh2.variables['mask_seq_sf_layer'][t0:t1,:] = masks.seq_sf_layer
        fh2.variables['mask_seq_sf_count'][t0:t1,:] = masks.seq_sf_count
        fh2.variables['mask_seq_dc_num'][t0:t1,:]   = masks.seq_dc_num
        fh2.variables['mask_seq_dc_layer'][t0:t1,:] = masks.seq_dc_layer
        fh2.variables['mask_seq_dc_count'][t0:t1,:] = masks.seq_dc_count
        fh2.variables['mask_gl_sf_layer'][t0:t1,:]  = masks.gl_sf_layer
        fh2.variables['mask_ml_sf_layer'][t0:t1,:]  = masks.ml_sf_layer
        fh2.variables['mask_gl_sf'][t0:t1,:]        = masks.gl_sf
        fh2.variables['mask_ml_sf'][t0:t1,:]        = masks.ml_sf
        fh2.variables['mask_gl_dc_layer'][t0:t1,:]  = masks.gl_dc_layer
        fh2.variables['mask_ml_dc_layer'][t0:t1,:]  = masks.ml_dc_layer
        fh2.variables['mask_gl_dc'][t0:t1,:]        = masks.gl_dc
        fh2.variables['mask_ml_dc'][t0:t1,:]        = masks.ml_dc
   
        # mixed layer characteristics
        fh2.variables['ml_h'][t0:t1,:]              = ml.height
        fh2.variables['ml_p'][t0:t1,:]              = ml.p
        fh2.variables['ml_T'][t0:t1,:]              = ml.T
        fh2.variables['ml_S'][t0:t1,:]              = ml.S
        fh2.variables['ml_Tu'][t0:t1,:]             = ml.Tu
        fh2.variables['ml_R'][t0:t1,:]              = ml.R
        fh2.variables['ml_r'][t0:t1,:]              = ml.r
        fh2.variables['ml_dS'][t0:t1,:]             = ml.dS
        fh2.variables['ml_dT'][t0:t1,:]             = ml.dT
        fh2.variables['ml_dr'][t0:t1,:]             = ml.dr
        fh2.variables['ml_h'][t0:t1,:]              = ml.height
        fh2.variables['ml_count'][t0:t1,:]          = ml.count
          
        #gradient layer characteristics
        fh2.variables['gl_dT'][t0:t1,:]             = gl.dT
        fh2.variables['gl_dS'][t0:t1,:]             = gl.dS
        fh2.variables['gl_dTdz'][t0:t1,:]           = gl.dTdz
        fh2.variables['gl_dSdz'][t0:t1,:]           = gl.dSdz
        fh2.variables['gl_dr'][t0:t1,:]             = gl.dr
        fh2.variables['gl_h'][t0:t1,:]              = gl.dist
        fh2.variables['gl_Tu'][t0:t1,:]             = gl.Tu
        fh2.variables['gl_R'][t0:t1,:]              = gl.R
        fh2.variables['gl_count'][t0:t1,:]          = gl.count
        fh2.close()
         
