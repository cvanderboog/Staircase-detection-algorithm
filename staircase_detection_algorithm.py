"""
Detection algorithm to find thermohaline staircases in vertical temperature and salinity profiles
A detailed description of the algorithm can be found in:

van der Boog, C.G. et al. (20xx), Global dataset of thermohaline staircases obtained from Argo
floats and Ice Tethered Profilers. Submitted to Earth System Science Data

Carine van der Boog (2020)
"""

import gsw
import numpy as np
import scipy.ndimage as ndimage
from scipy.signal import argrelextrema

def get_mixed_layers(p,ct,sa,c1,c2,c3,c4):
  copyright()
  """
  all data should be on 1 dbar resolution
  p = pressure (dbar)
  ct = conservative temperature (degrees celsius)
  sa = absolute salinity (g kg-1)

  c1 = mixed layer gradient with gradient detection method
  c2 = density difference of mixed layer
  c3 = averaging window to obtain background profiles
  c4 = ratio between mixed layer height and gradient layer height

  """
  # maximum mixed layers per profiles = 100
  max_count = 200 

  """
  step 0:  define classes
  """
  class ml:
    pass

  class gl:
    pass

  class masks:
    pass

  empty_ml_array = np.ma.array(np.ma.zeros((len(p[:,0]),max_count)),mask=np.ones((len(p[:,0]),max_count)))
  empty_gl_array = np.ma.array(np.ma.zeros((len(p[:,0]),max_count-1)),mask=np.ones((len(p[:,0]),max_count-1)))

  """
  step 1a:  Mixed layers - compute mixed layers
  """
  # compute vertical gradients 
  rho            = gsw.sigma1(sa,ct)
  dTdz           = cent_derivative2d(ct,p)
  dSdz           = cent_derivative2d(sa,p)
  drdz           = cent_derivative2d(rho,p)
  
  # compute density gradients    - it is assumed that the alpha  
  T_av           = moving_average2d(ct,c3)      #c3
  S_av           = moving_average2d(sa,c3)      #c3
  alpha          = gsw.alpha(S_av,T_av,p)
  beta           = gsw.beta(S_av,T_av,p)
  dTdz_av        = cent_derivative2d(T_av,p)
  dSdz_av        = cent_derivative2d(S_av,p)

  T_av2          = moving_average2d(ct,50)      #c3
  S_av2          = moving_average2d(sa,50)      #c3
  dTdz_av2       = cent_derivative2d(T_av2,p)
  dSdz_av2       = cent_derivative2d(S_av2,p)
  Tu             = np.arctan2(-alpha*dTdz_av2-beta*dSdz_av2,-alpha*dTdz_av2+beta*dSdz_av2)*360/(2*np.pi)
  R              = alpha*dTdz_av2/(beta*dSdz_av2)

  # mask where there is no moving average
  dTdz.mask      = dTdz_av.mask
  dSdz.mask      = dSdz_av.mask
  sa.mask        = S_av.mask
  ct.mask        = T_av.mask

  # make mask where only mixed layers are detected
  gradientT      = np.ma.masked_outside(-alpha*1028*dTdz,-c1,c1)    #c1
  gradientS      = np.ma.masked_outside(beta*1028*dSdz,-c1,c1)      #c1
  gradientrho    = np.ma.masked_outside(drdz,-c1,c1)                #c1
  layer_mask     = gradientT.mask + gradientS.mask + gradientrho.mask 

  """
  step 1b:  get mixed layer extent
  """
  # count mixed layers           - start new mixed layer by 2 not counts
  p_new          = np.ma.zeros(p.shape,dtype=np.int32)
  p_new.mask     = np.ma.copy(layer_mask)
  for k in range(len(p_new[:,0])):
    start = 1
    if p[k,:].max():
      p_new[k,0] = 0
      for l in range(len(p_new[0,:])-1):
        if p_new.mask[k,l]==False:
          p_new[k,l] = start
        if p_new.mask[k,l]==False and p_new.mask[k,l+1]==True and p_new.mask[k,l+1]==True:
         start   = start+1
  ml_count_tmp       = p_new.max(axis=1)

  # loop through each mixed layer
  for k in range(len(ml_count_tmp)):
    if ml_count_tmp[k]:
      for l in range(np.min([max_count,ml_count_tmp[k]])):
        # locate the layer
        layer = np.where(p_new[k,:]==l+1)[0][:]
        ml_T      = np.ma.mean(-alpha[k*np.ones(len(layer),dtype=np.int32),layer]*1028*ct[k*np.ones(len(layer),dtype=np.int32),layer])
        ml_S      = np.ma.mean(beta[k*np.ones(len(layer),dtype=np.int32),layer]*1028*sa[k*np.ones(len(layer),dtype=np.int32),layer])
        ml_r      = np.ma.mean(rho[k*np.ones(len(layer),dtype=np.int32),layer])
        temp_prof      = np.ma.copy(-alpha[k,:]*1028*ct[k,:])
        salt_prof      = np.ma.copy(beta[k,:]*1028*sa[k,:])
        rho_prof       = np.ma.copy(rho[k,:])

        # find the layers within the temperature and salinity range
        ml_temp        = np.ma.masked_outside(temp_prof,ml_T-0.5*c2,ml_T+0.5*c2)    #c2
        ml_salt        = np.ma.masked_outside(salt_prof,ml_S-0.5*c2,ml_S+0.5*c2)    #c2
        ml_rho         = np.ma.masked_outside(rho_prof,ml_r-0.5*c2,ml_r+0.5*c2)    #c2
        mixed_layers   = np.where(ml_temp.mask+ml_salt.mask+ml_rho.mask==False)[0][:]
        new_indices    = np.sort(np.unique(np.concatenate((layer,mixed_layers),axis=0)))

        #check if new indices are continuous
        if len(new_indices)>1:
          if np.max(new_indices[1:]-new_indices[:-1])>1:
            split          = np.where(np.diff(new_indices) != 1)[0][:]+1
            if len(split)==1:
              if split < 0.5*len(new_indices):
                new_indices = new_indices[split[0]:]
              else:
                new_indices = new_indices[:split[0]]
            else:
              #find longest sequence
              new           = np.argmax(np.append(np.append(split[0],split[1:]-split[:-1]),len(new_indices)-split[-1]))
              if new >= len(split):
                new_indices   = new_indices[split[-1]:]
              elif new == 0:
                new_indices   = new_indices[:split[new]]
              else:
                new_indices   = new_indices[split[new-1]:split[new]]
        layer_mask[k*np.ones(len(new_indices),dtype=np.int32),new_indices] = False   #update mask
        layer_mask[k,new_indices[-1]+1] = True            #mask top and bottom of layer                #MASK GL
        layer_mask[k,new_indices[0]-1] = True            #mask top and bottom of layer                 #MASK GL

  # update the mixed layer count 
  p_new       = np.ma.zeros(p.shape,dtype=np.int32)
  p_new.mask  = np.ma.copy(layer_mask)

  for k in range(len(p_new[:,0])):
    start = 1
    if p[k,:].max():
      p_new[k,0] = 0
      for l in range(len(p_new[0,:])-1):
        if p_new.mask[k,l]==False:
          p_new[k,l] = start
        if p_new.mask[k,l]==False and p_new.mask[k,l+1]==True:
          start = start+1
  ml_count_tmp    = p_new.max(axis=1) 

  """
  step 1c: mixed layer and gradient layer properties
  """
  # gradient layers are defined as layer between two mixed layers. 
  # variables that need to be computed per mixed layer
  ml.T       = np.ma.copy(empty_ml_array)   # average temperature per mixed layer
  ml.S       = np.ma.copy(empty_ml_array)   # average salinity mixed layer
  ml.r       = np.ma.copy(empty_ml_array)   # average density mixed layer
  ml.p       = np.ma.copy(empty_ml_array)   # average pressure mixed layer
  ml.dT      = np.ma.copy(empty_ml_array)   # temperature difference within mixed layer
  ml.dS      = np.ma.copy(empty_ml_array)   # salinity difference within mixed layer
  ml.dr      = np.ma.copy(empty_ml_array)   # density difference within mixed layer
  ml.height  = np.ma.copy(empty_ml_array)   # mixed layer height
  ml.Tu      = np.ma.copy(empty_ml_array)   # average Turner angle mixed layer
  ml.R       = np.ma.copy(empty_ml_array)   # average density ratio of the mixed layer
  ml.count   = np.ma.array(np.ma.copy(empty_ml_array),dtype=np.int32)
                                            # mixed layer number (counted from top)
  pmin     = np.ma.array(np.ma.zeros((len(ml_count_tmp),max_count)),mask=np.ones((len(ml_count_tmp),max_count)),dtype=np.int32)
  pmax     = np.ma.array(np.ma.zeros((len(ml_count_tmp),max_count)),mask=np.ones((len(ml_count_tmp),max_count)),dtype=np.int32)

  gl.dTdz  = np.ma.copy(empty_gl_array)     # temperature difference within gradient layer
  gl.dSdz  = np.ma.copy(empty_gl_array)     # temperature difference within gradient layer
  gl.dT    = np.ma.copy(empty_gl_array)     # temperature difference within gradient layer
  gl.dS    = np.ma.copy(empty_gl_array)     # salinity difference within gradient layer
  gl.dr    = np.ma.copy(empty_gl_array)     # density difference within gradient layer
  gl.dist  = np.ma.copy(empty_gl_array)     # distance between gradient layer
  gl.Tu    = np.ma.copy(empty_gl_array)     # average turner angle gradient layuer
  gl.R     = np.ma.copy(empty_gl_array)     # average density ratio gradient layer
  gl.count = np.ma.array(np.ma.copy(empty_gl_array),dtype=np.int32)
                                            # gradient layer number (counted from top)
  gl.intersect = np.ma.array(np.ma.copy(empty_gl_array),dtype=np.int32) # slope changes
  # masks for ct and sa where the layers are  
  masks.gl_num  = np.ma.array(np.zeros(p.shape),mask=np.ones(p.shape),dtype=np.int32)
  masks.ml_num  = np.ma.array(np.zeros(p.shape),mask=np.ones(p.shape),dtype=np.int32)

  # loop through each mixed layer
  for k in range(len(ml_count_tmp)):
    if ml_count_tmp[k]:
      for l in range(np.min([max_count,ml_count_tmp[k]])):
        layer = np.where(p_new[k,:]==l+1)[0][:]
        ks             = k*np.ones(len(layer),dtype=np.int32)
        ml.height[k,l] = len(np.arange(layer.min(),layer.max()+1))
        ml.p[k,l]      = p[ks,layer].mean()
        ml.T[k,l]      = ct[ks,layer].mean()
        ml.S[k,l]      = sa[ks,layer].mean()
        ml.r[k,l]      = rho[ks,layer].mean()
        ml.R[k,l]      = R[ks,layer].mean()
        ml.Tu[k,l]      = Tu[ks,layer].mean()
        ml.dr[k,l]     = rho[ks,layer].max()-rho[ks,layer].min()
        rhomax         = rho[ks,layer].argmax()
        rhomin         = rho[ks,layer].argmin()
        ml.dT[k,l]     = ct[k,layer[rhomax]]-ct[k,layer[rhomin]]
        ml.dS[k,l]     = sa[k,layer[rhomax]]-sa[k,layer[rhomin]]

        ml.count[k,l]     = l+1
        masks.ml_num[ks,layer] = l+1

        pmin[k,l]       = np.int(layer.min())
        pmax[k,l]       = np.int(layer.max())

      for g in range(np.min([max_count,ml_count_tmp[k]])-1):
        indices        = np.arange(pmax[k,g]+1,pmin[k,g+1])
        ks             = k*np.ones(len(indices),dtype=np.int32)
        gl.dT[k,g]     = ct[ks,indices].max()-ct[ks,indices].min()
        gl.dS[k,g]     = sa[ks,indices].max()-sa[ks,indices].min()
        gl.dr[k,g]     = rho[ks,indices].max()-rho[ks,indices].min()
        gl.dist[k,g]   = len(indices)
        gl.Tu[k,g]     = Tu[ks,indices].mean()
        gl.R[k,g]      = R[ks,indices].mean()
        gl.count[k,g]  = g+1
        gl.intersect[k,g] = len(argrelextrema(ct[ks,indices], np.greater)[0])+len(argrelextrema(ct[ks,indices], np.less)[0]) 
        gl.dTdz[k,g]   = dTdz_av2[ks,indices].mean()
        gl.dSdz[k,g]   = dSdz_av2[ks,indices].mean()
        masks.gl_num.mask[ks,indices] = False
        masks.gl_num[ks,indices] = g+1

  """
  step 2+3:  label each gradient layer that fullfills the requirements for the gradient layer height
  and gradient layer properties (temperature, salinity and density variations) 
  """
  ml_h_max  = np.ma.array([np.ma.abs(ml.height[:,1:]),np.ma.abs(ml.height[:,:-1])]).min(axis=0)
  # check whether this makes a difference
  ml_h_max2 = np.ma.copy(ml_h_max)
  ml_h_max2[ml_h_max2>c4] = c4

  ml_dS_max = np.ma.array([np.ma.abs(ml.dS[:,1:]),np.ma.abs(ml.dS[:,:-1])]).max(axis=0)
  ml_dT_max = np.ma.array([np.ma.abs(ml.dT[:,1:]),np.ma.abs(ml.dT[:,:-1])]).max(axis=0)
  ml_dr_max = np.ma.array([np.ma.abs(ml.dr[:,1:]),np.ma.abs(ml.dr[:,:-1])]).max(axis=0)

  masks.qc = np.ma.array(0*np.ma.copy(gl.count)+1,dtype=np.int32)
  masks.qc[np.ma.abs(ml_dT_max) < np.ma.abs(gl.dT)]   +=2
  masks.qc[np.ma.abs(ml_dS_max) < np.ma.abs(gl.dS)]   +=4
  masks.qc[np.ma.abs(ml_dr_max) < np.ma.abs(gl.dr)]   +=8
  masks.qc[np.ma.abs(ml_h_max)  > np.ma.abs(gl.dist)] +=16
  masks.qc[np.ma.abs(ml_h_max2) > np.ma.abs(gl.dist)] +=32
  masks.qc[gl.intersect < 3]                          +=64
 

  """
  step 4:  flag each gradient layers as 'salt-finger' or a 'diffusive convective' or nothing
  """
  masks.regime = np.ma.array(np.ma.copy(empty_gl_array),dtype=np.int32)

  for k in range(len(ml_count_tmp)):
    ml_diffT = ml.T[k,1:]-ml.T[k,:-1]
    ml_diffS = ml.S[k,1:]-ml.S[k,:-1]

    ml_diffT[(ml_diffT<0) & (ml_diffT.mask==False)] = -1
    ml_diffS[(ml_diffS<0) & (ml_diffS.mask==False)] = -1
    ml_diffT[(ml_diffT>0) & (ml_diffT.mask==False)] = 1
    ml_diffS[(ml_diffS>0) & (ml_diffS.mask==False)] = 1

    ml_diff = ml_diffT+ml_diffS
    masks.regime[k,:] = ml_diff

  """
  step 5:  select sequences of mixed layers
  """
  #salt finger regime
  gl_tmp_sf                        = np.ma.array(np.copy(masks.qc),mask=np.copy(masks.qc.mask))
  gl_tmp_sf.mask[gl_tmp_sf<127]    = True
  gl_tmp_sf.mask[masks.regime!=-2] = True

  masks.seq_sf_num              = np.ma.array(np.ma.copy(empty_gl_array),dtype=np.int32)
  masks.seq_sf_layer            = np.ma.array(np.ma.copy(empty_gl_array),dtype=np.int32)
  masks.seq_sf_count            = np.ma.array(np.ma.copy(empty_gl_array),dtype=np.int32)

  for k in range(len(ml_count_tmp)):
    seq_count_tmp_sf = 1
    gl_count_tmp_sf  = 1
    for l in range(max_count-1):
      if gl_tmp_sf.mask[k,l]==False:
        masks.seq_sf_num[k,l] = seq_count_tmp_sf
        masks.seq_sf_layer[k,l] = gl_count_tmp_sf
        gl_count_tmp_sf = gl_count_tmp_sf+1
      if l > 0:
        if gl_tmp_sf.mask[k,l] == True and gl_tmp_sf.mask[k,l-1] == False:
          indices = np.where(masks.seq_sf_num[k,:]==seq_count_tmp_sf)[0][:]
          masks.seq_sf_count[k*np.ones(len(indices),dtype=np.int32),indices]=gl_count_tmp_sf-1
          seq_count_tmp_sf = seq_count_tmp_sf+1
          gl_count_tmp_sf  = 1

  #diffusive convective regime
  gl_tmp_dc                       = np.ma.array(np.copy(masks.qc),mask=np.copy(masks.qc.mask))
  gl_tmp_dc.mask[gl_tmp_dc<127]   = True
  gl_tmp_dc.mask[masks.regime!=2] = True

  masks.seq_dc_num              = np.ma.array(np.ma.copy(empty_gl_array),dtype=np.int32)
  masks.seq_dc_layer            = np.ma.array(np.ma.copy(empty_gl_array),dtype=np.int32)
  masks.seq_dc_count            = np.ma.array(np.ma.copy(empty_gl_array),dtype=np.int32)

  for k in range(len(ml_count_tmp)):
    seq_count_tmp_dc = 1
    gl_count_tmp_dc  = 1
    for l in range(max_count-1):
      if gl_tmp_dc.mask[k,l]==False:
        masks.seq_dc_num[k,l] = seq_count_tmp_dc
        masks.seq_dc_layer[k,l] = gl_count_tmp_dc
        gl_count_tmp_dc = gl_count_tmp_dc+1
      if l > 0:
        if gl_tmp_dc.mask[k,l] == True and gl_tmp_dc.mask[k,l-1] == False:
          indices = np.where(masks.seq_dc_num[k,:]==seq_count_tmp_dc)[0][:]
          masks.seq_dc_count[k*np.ones(len(indices),dtype=np.int32),indices]=gl_count_tmp_dc-1
          seq_count_tmp_dc = seq_count_tmp_dc+1
          gl_count_tmp_dc  = 1

  """
  step 9:  mark the sequences in an array
  """
  masks.gl_sf_layer = np.ma.array(np.ma.copy(empty_gl_array),dtype=np.int32)
  masks.ml_sf_layer = np.ma.array(np.ma.copy(empty_ml_array),dtype=np.int32)
  masks.gl_sf       = np.ma.array(np.zeros(p.shape),mask=np.ones(p.shape),dtype=np.int32)
  masks.ml_sf       = np.ma.array(np.zeros(p.shape),mask=np.ones(p.shape),dtype=np.int32)
  masks.gl_dc_layer = np.ma.array(np.ma.copy(empty_gl_array),dtype=np.int32)
  masks.ml_dc_layer = np.ma.array(np.ma.copy(empty_ml_array),dtype=np.int32)
  masks.gl_dc       = np.ma.array(np.zeros(p.shape),mask=np.ones(p.shape),dtype=np.int32)
  masks.ml_dc       = np.ma.array(np.zeros(p.shape),mask=np.ones(p.shape),dtype=np.int32)

  for k in range(len(ml_count_tmp)):
    # salt finger regime
    if masks.seq_sf_num[k,:].max():
      for seq in range(masks.seq_sf_num[k,:].max()):
        indx_seq_gl = np.where(masks.seq_sf_num[k,:]==seq+1)[0][:]
        indx_seq_ml = np.append(indx_seq_gl,indx_seq_gl[-1]+1)

        ks_gl = k*np.ones(len(indx_seq_gl),dtype=np.int32)
        mask_gl_tmp = np.ma.copy(masks.gl_num[k,:])
        indx_seq_gl_long = np.where(np.ma.masked_outside(mask_gl_tmp,gl.count[ks_gl,indx_seq_gl].min(),gl.count[ks_gl,indx_seq_gl].max()).mask==False)[0][:]
        ks_gl_long =  k*np.ones(len(indx_seq_gl_long),dtype=np.int32)

        ks_ml = k*np.ones(len(indx_seq_ml),dtype=np.int32)
        mask_ml_tmp = np.ma.copy(masks.ml_num[k,:])
        indx_seq_ml_long = np.where(np.ma.masked_outside(mask_ml_tmp,ml.count[ks_ml,indx_seq_ml].min(),ml.count[ks_ml,indx_seq_ml].max()).mask==False)[0][:]
        ks_ml_long =  k*np.ones(len(indx_seq_ml_long),dtype=np.int32)

        masks.gl_sf_layer[ks_gl,indx_seq_gl]           = len(indx_seq_gl)
        masks.ml_sf_layer[ks_ml,indx_seq_ml]           = len(indx_seq_ml)
        masks.gl_sf[ks_gl_long,indx_seq_gl_long]       = len(indx_seq_gl)
        masks.ml_sf[ks_ml_long,indx_seq_ml_long]       = len(indx_seq_ml) 

    # diffusive convection regime
    if masks.seq_dc_num[k,:].max():
      for seq in range(masks.seq_dc_num[k,:].max()):
        indx_seq_gl = np.where(masks.seq_dc_num[k,:]==seq+1)[0][:]
        indx_seq_ml = np.append(indx_seq_gl,indx_seq_gl[-1]+1)

        ks_gl = k*np.ones(len(indx_seq_gl),dtype=np.int32)
        mask_gl_tmp = np.ma.copy(masks.gl_num[k,:])
        indx_seq_gl_long = np.where(np.ma.masked_outside(mask_gl_tmp,gl.count[ks_gl,indx_seq_gl].min(),gl.count[ks_gl,indx_seq_gl].max()).mask==False)[0][:]
        ks_gl_long =  k*np.ones(len(indx_seq_gl_long),dtype=np.int32)

        ks_ml = k*np.ones(len(indx_seq_ml),dtype=np.int32)
        mask_ml_tmp = np.ma.copy(masks.ml_num[k,:])
        indx_seq_ml_long = np.where(np.ma.masked_outside(mask_ml_tmp,ml.count[ks_ml,indx_seq_ml].min(),ml.count[ks_ml,indx_seq_ml].max()).mask==False)[0][:]
        ks_ml_long =  k*np.ones(len(indx_seq_ml_long),dtype=np.int32)

        masks.gl_dc_layer[ks_gl,indx_seq_gl]           = len(indx_seq_gl)
        masks.ml_dc_layer[ks_ml,indx_seq_ml]           = len(indx_seq_ml) 
        masks.gl_dc[ks_gl_long,indx_seq_gl_long]       = len(indx_seq_gl)
        masks.ml_dc[ks_ml_long,indx_seq_ml_long]       = len(indx_seq_ml) 
  return gl, ml, masks



def cent_derivative2d(f,z):
  dz1     = z[:,1:-1]-z[:,:-2]
  dz2     = z[:,2:]-z[:,1:-1]
  fkp1    = f[:,2:]
  fkm1    = f[:,:-2]
  fk      = f[:,1:-1]

  dfdz    = np.ma.zeros(f.shape)
  dfdz[:,1:-1] = ( dz1**2*fkp1 + ( dz2**2-dz1**2 )*fk - dz2**2*fkm1
                 ) / ( dz1*dz2*(dz1+dz2) )
  dfdz.mask[dfdz==0]=True 
  return dfdz

def moving_average2d(dataset,window):
  start = np.argmin(np.abs(np.cumsum(np.ma.logical_not(dataset.mask),axis=1)-1),axis=1)
  end   = np.argmax(np.ma.array(np.ma.cumsum(np.ma.logical_not(dataset.mask),axis=1),mask=dataset.mask),axis=1)
  index = np.int(window/2)+1
  for k in range(len(start)):
    dataset[k,:index+start[k]+1].mask = True
    dataset[k,end[k]-index-1:].mask = True
  weights = np.ma.repeat(1.0,window)/window
  mav = np.ma.array(ndimage.convolve1d(dataset,weights,axis=1,mode='constant',cval=0),mask=dataset.mask)
  return mav

def copyright():
  print('___________________')
  print('Algorithm to detect thermohaline staircases from Argo floats and Ice tethered profilers')
  print('Copyright (C) 2020  Carine van der Boog')
  print('')

  print('This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation version 3 of the License. \n This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. \n You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.')
  print('')
  print('___________________')
  print('')
  print('')
  print('')
  print('')
  return
