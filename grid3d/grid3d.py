#!/usr/bin/env python3
# grid3d.py

import os
import numpy as np
import h5py
#import csv
import picdefs

class Grid3d:
  def __init__(self, piccode):
    self.piccode = piccode
  
  def read(self, file):
    self.file = file
    print('Reading: ', self.file)
    if self.piccode == picdefs.code.HIPACE:
      self.read_hipace()
    elif self.piccode == picdefs.code.OSIRIS:
      self.read_osiris()
    print('Read-in completed.')
  def read_hipace(self): 

    fpath, fext = os.path.splitext(self.file)

    if any(fext == s for s in picdefs.fexts.hdf5):
      self.read_hipace_hdf5()
    else:
      print('ERROR:\tExtension of file "%s" is not supported!' %fname)
      print('\tAllowed file extensions: ',list(picdefs.fexts.hdf5))
      sys.exit()

  def read_hipace_hdf5(self):
    with h5py.File(self.file,'r') as hf:
      
      # Reading dataset
      self.h5keys = list(hf.keys())
      if len(self.h5keys) == 0:
        print('ERROR:\tHDF5 file "%s" does not contain any dataset!' %self.file)
        sys.exit()
      elif len(self.h5keys) == 1:
        self.data = np.array(hf.get(self.h5keys[0]))
      else:
        print('ERROR:\tHDF5 file "%s" contains more than one dataset!' %self.file)
        sys.exit()

      # Printing attributes
      for item in hf.attrs.keys():
        print(item + ":", hf.attrs[item])
      
      # Reading attributes 
      self.nx = hf.attrs[   picdefs.hipace.h5attrkeys.nx ] 
      self.xmin = hf.attrs[ picdefs.hipace.h5attrkeys.xmin ]
      self.xmax = hf.attrs[ picdefs.hipace.h5attrkeys.xmax ]
      self.time = hf.attrs[ picdefs.hipace.h5attrkeys.time ]
      self.dt = hf.attrs[   picdefs.hipace.h5attrkeys.dt ]
      self.type = hf.attrs[ picdefs.hipace.h5attrkeys.type ]
      self.name = hf.attrs[ picdefs.hipace.h5attrkeys.name ]
      
  def get_x_arr(self,dim):
    return np.linspace(self.xmin[dim],self.xmax[dim],self.nx[dim])
  def get_z_arr(self):  
    return np.linspace(self.time+self.xmin[0],self.time+self.xmax[0],self.nx[0])
  def get_zeta_arr(self):  
    return np.linspace(self.xmin[0],self.xmax[0],self.nx[0])
  def get_xi_arr(self):  
    return np.linspace(-self.xmin[0],-self.xmax[0],self.nx[0])  
  def get_nt(self):
    return round(self.time/self.dt)
   
  def read_osiris(self):
    print('OSIRIS Read-in not yet implemented!')
    sys.exit()
    