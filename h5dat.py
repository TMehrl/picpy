#!/usr/bin/env python3
# h5dat.py

import os
import numpy as np
import h5py
import picdefs
import sys


def hdf5_check(file):
  fpath, fext = os.path.splitext(file)
  if any(fext == s for s in picdefs.fexts.hdf5):
    pass
  else:
    print('Error:\tExtension of file "%s" is not supported!' %fname)
    print('\tAllowed file extensions: ',list(picdefs.fexts.hdf5))
    sys.exit()   

def check_if_h5file_is_read(if_read):
  if if_read: 
    pass
  else:
    print('Error: HDF5 file not yet read!')
    sys.exit() 


def print_datasets(file):
  with h5py.File(file,'r') as hf:
    # Printing attributes
    print('HDF5 file datasets:')
    for item in hf.keys():
      print('\t' + item + ":", hf[item])

def print_attributes(file):
  with h5py.File(file,'r') as hf:
    # Printing attributes
    print('HDF5 file attributes:')
    for item in hf.attrs.keys():
      print('\t' + item + ":", hf.attrs[item])


class RAW:
  def __init__(self, piccode):
    self.piccode = piccode
    self.if_h5file_read = False

  def print_attributes(self, file):
    print_attributes(file)

  def print_datasets(self, file):
    print_datasets(file)

  def read(self, file):
    hdf5_check(file) 
    self.file = file  
    if self.piccode == picdefs.code.hipace:
      self.read_hipace()
    elif self.piccode == picdefs.code.osiris:
      self.read_osiris()

  def read_hipace(self): 
  
    with h5py.File(self.file,'r') as hf:
      
      # Reading datasets
      self.x1 = np.array(hf.get(  picdefs.hipace.h5.rawkeys.x1 ))
      self.x2 = np.array(hf.get(  picdefs.hipace.h5.rawkeys.x2 ))
      self.x3 = np.array(hf.get(  picdefs.hipace.h5.rawkeys.x3 ))
      self.q =  np.array(hf.get(  picdefs.hipace.h5.rawkeys.q  ))
      self.p1 = np.array(hf.get(  picdefs.hipace.h5.rawkeys.p1 ))    
      self.p2 = np.array(hf.get(  picdefs.hipace.h5.rawkeys.p2 ))
      self.p3 = np.array(hf.get(  picdefs.hipace.h5.rawkeys.p3 ))
            
      # Reading attributes 
      self.nx = hf.attrs[   picdefs.hipace.h5.attrkeys.nx ] 
      self.xmin = hf.attrs[ picdefs.hipace.h5.attrkeys.xmin ]
      self.xmax = hf.attrs[ picdefs.hipace.h5.attrkeys.xmax ]
      self.time = hf.attrs[ picdefs.hipace.h5.attrkeys.time ]
      self.dt = hf.attrs[   picdefs.hipace.h5.attrkeys.dt ]
      type_bytes = hf.attrs[ picdefs.hipace.h5.attrkeys.type ]
      name_bytes = hf.attrs[ picdefs.hipace.h5.attrkeys.name ]
      self.type = type_bytes[0].decode('UTF-8')
      self.name = name_bytes[0].decode('UTF-8')
    
    self.if_h5file_read = True

  def read_osiris(self):
    print('Error: OSIRIS Read-in not yet implemented!')
    sys.exit()

#### 3D-grid 
class Grid3d:
  def __init__(self, piccode):
    self.piccode = piccode
    self.if_h5file_read = False

  def print_attributes(self, file):
    print_attributes(file)

  def print_datasets(self, file):
    print_datasets(file)

  def read(self, file):
    hdf5_check(file) 
    self.file = file  
    if self.piccode == picdefs.code.hipace:
      self.read_hipace()
    elif self.piccode == picdefs.code.osiris:
      self.read_osiris()

  
  def read_osiris(self):
    print('Error: OSIRIS Read-in not yet implemented!')
    sys.exit()
      
  def read_hipace(self): 
  
    with h5py.File(self.file,'r') as hf:
      
      # Reading dataset (here not caring how dataset is called)
      self.h5keys = list(hf.keys())
      if len(self.h5keys) == 0:
        print('Error:\tHDF5 file "%s" does not contain any dataset!' %(self.file) )
        sys.exit()
      elif len(self.h5keys) == 1:
        self.data = np.array(hf.get(self.h5keys[0]))
      else:
        print('Error:\tHDF5 file "%s" contains more than one dataset!' %(self.file) )
        sys.exit()
      
      # Reading attributes 
      self.nx = hf.attrs[   picdefs.hipace.h5.attrkeys.nx ] 
      self.xmin = hf.attrs[ picdefs.hipace.h5.attrkeys.xmin ]
      self.xmax = hf.attrs[ picdefs.hipace.h5.attrkeys.xmax ]
      self.time = hf.attrs[ picdefs.hipace.h5.attrkeys.time ]
      self.dt = hf.attrs[   picdefs.hipace.h5.attrkeys.dt ]
      type_bytes = hf.attrs[ picdefs.hipace.h5.attrkeys.type ]
      name_bytes = hf.attrs[ picdefs.hipace.h5.attrkeys.name ]
      self.type = type_bytes[0].decode('UTF-8')
      self.name = name_bytes[0].decode('UTF-8')
    
    self.if_h5file_read = True
        
  def get_x_arr(self,dim):
    check_if_h5file_is_read(self.if_h5file_read)
    if self.piccode == picdefs.code.hipace:
      return np.linspace(self.xmin[dim],self.xmax[dim],self.nx[dim])
    else:
      print('Error: OSIRIS part not yet implemented!')
      sys.exit()
      
  def get_z_arr(self):
    check_if_h5file_is_read(self.if_h5file_read) 
    if self.piccode == picdefs.code.hipace:
      return np.linspace(self.time+self.xmin[0],self.time+self.xmax[0],self.nx[0])
    else:
      print('Error: OSIRIS part not yet implemented!')
      sys.exit()  
            
  def get_zeta_arr(self):
    check_if_h5file_is_read(self.if_h5file_read) 
    if self.piccode == picdefs.code.hipace:
      return np.linspace(self.xmin[0],self.xmax[0],self.nx[0])
    else:
      print('Error: OSIRIS part not yet implemented!')
      sys.exit()
            
  def get_xi_arr(self):
    check_if_h5file_is_read(self.if_h5file_read) 
    if self.piccode == picdefs.code.hipace:
      return np.linspace(-self.xmin[0],-self.xmax[0],self.nx[0])
    else:
      print('Error: OSIRIS part not yet implemented!')
      sys.exit()
                  
  def get_nt(self):
    check_if_h5file_is_read(self.if_h5file_read)       
    if self.piccode == picdefs.code.hipace:
      return round(self.time/self.dt)
    else:
      print('Error: OSIRIS part not yet implemented!')
      sys.exit()
    