#!/usr/bin/env python3
# h5dat.py

import os
import numpy as np
import h5py
import picdefs
import sys


class H5File:
  def __init__(self, file, check=False):

    # Allowed hdf5 extensions:
    self.__h5exts = ['.h5','.hdf5']

    self.file = file
    
    if check:
      self.hdf5_check()

  def hdf5_check(self):
    fpath, fext = os.path.splitext(self.file)
    if any(fext == s for s in self.__h5exts):
      pass
    else:
      print('Error:\tExtension of file "%s" is not supported!' %fname)
      print('\tAllowed file extensions: ',list(self.__h5exts))
      sys.exit()

  def if_hdf5_file(self):
    fpath, fext = os.path.splitext(self.file)
    if any(fext == s for s in self.__h5exts):
      return True
    else:
      return False

  def print_datasets(self):
    with h5py.File(self.file,'r') as hf:
      # Printing attributes
      print('HDF5 file datasets:')
      for item in hf.keys():
        print('\t' + item + ":", hf[item])

  def print_attributes(self):
    with h5py.File(self.file,'r') as hf:
      # Printing attributes
      print('HDF5 file attributes:')
      for item in hf.attrs.keys():
        print('\t' + item + ":", hf.attrs[item])

  def get_allowed_h5exts(self):
    return self.__h5exts


class H5PIC(H5File):
  def __init__(self, file, piccode):

    H5File.__init__(self, file, check=True)

    self.__hipace = 'hipace'
    self.__osiris = 'osiris'

    self.piccode = piccode

    if self.piccode == self.__hipace:
      self.read_hipace()
    elif self.piccode == self.__osiris:
      self.read_osiris()

    
class RAW(H5PIC):
  def __init__(self, file, piccode=picdefs.code.hipace):
    H5PIC.__init__(self, file, piccode)
    
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
      
      self.npart = len(self.q)
            
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

  def read_osiris(self):
    print('Error: OSIRIS Read-in not yet implemented!')
    sys.exit()


#### 3D-grid 
class Grid3d(H5PIC):
  def __init__(self, file, piccode=picdefs.code.hipace):
    H5PIC.__init__(self, file, piccode)

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
        
  def get_x_arr(self,dim):
    if self.piccode == picdefs.code.hipace:
      return np.linspace(self.xmin[dim],self.xmax[dim],self.nx[dim])
    else:
      print('Error: OSIRIS part not yet implemented!')
      sys.exit()
      
  def get_z_arr(self):
    if self.piccode == picdefs.code.hipace:
      return np.linspace(self.time+self.xmin[0],self.time+self.xmax[0],self.nx[0])
    else:
      print('Error: OSIRIS part not yet implemented!')
      sys.exit()  
            
  def get_zeta_arr(self):
    if self.piccode == picdefs.code.hipace:
      return np.linspace(self.xmin[0],self.xmax[0],self.nx[0])
    else:
      print('Error: OSIRIS part not yet implemented!')
      sys.exit()
            
  def get_xi_arr(self): 
    if self.piccode == picdefs.code.hipace:
      return np.linspace(-self.xmin[0],-self.xmax[0],self.nx[0])
    else:
      print('Error: OSIRIS part not yet implemented!')
      sys.exit()
                  
  def get_nt(self):      
    if self.piccode == picdefs.code.hipace:
      return round(self.time/self.dt)
    else:
      print('Error: OSIRIS part not yet implemented!')
      sys.exit()


class SliceMoms(H5File):
  def __init__(self, file):    
    
    H5File.__init__(self, file, check=True)
    
    self.read()
    
  def read(self):        
    with h5py.File(self.file,'r') as hf:
      # read order of moments
      order = 2
      # Reading datasets
      self.time_array = np.array(hf.get( 'time_array' ))
      self.zeta_array = np.array(hf.get( 'zeta_array' ))
      self.avgx1 = np.array(hf.get( 'avgx1' ))
      self.avgx2 = np.array(hf.get( 'avgx2' ))
      self.avgx3 = np.array(hf.get( 'avgx3' ))
      self.avgp1 = np.array(hf.get( 'avgp1' ))
      self.avgp2 = np.array(hf.get( 'avgp2' ))
      self.avgp3 = np.array(hf.get( 'avgp3' ))
      
      if order > 1:
        self.avgx1sq = np.array(hf.get( 'avgx1sq' ))
        self.avgx2sq = np.array(hf.get( 'avgx2sq' ))
        self.avgx3sq = np.array(hf.get( 'avgx3sq' ))
        self.avgp1sq = np.array(hf.get( 'avgp1sq' ))
        self.avgp2sq = np.array(hf.get( 'avgp2sq' ))
        self.avgp3sq = np.array(hf.get( 'avgp3sq' ))
        self.avgx1p1 = np.array(hf.get( 'avgx1p1' ))
        self.avgx2p2 = np.array(hf.get( 'avgx2p2' ))
        self.avgx3p3 = np.array(hf.get( 'avgx3p3' ))


class DIR:
  def __init__(self, dir):
    if os.path.isdir(dir):
      sys.stdout.write('Reading directory: %s\n' % dir)
      sys.stdout.flush()
      
      self.dir = dir
    else:
      print(dir + ' is not a directory!')
      sys.exit()
  
  def list_files(self, strpattern): 
    self.flist=[]
    ftime=[]
    flist_usort=[]
    for root, dirs, files in os.walk(self.dir):
      for name in files:
        h5file = H5File(name)
        if (strpattern in name) & h5file.if_hdf5_file():
          fpath, fext = os.path.splitext(name)
          ftime.append(int(fpath[-picdefs.hipace.h5.Nfname_digits:]))
          flist_usort.append(name)
      idx = np.argsort(np.asarray(ftime))
      self.nf = len(flist_usort)
      # Sorting
      for i in range(0,self.nf):
        self.flist.append(flist_usort[idx[i]])
        
  def filepath(self, i):
    fstr = '%s/%s' % (self.dir, self.flist[i])
    return fstr