#!/usr/bin/env python3
# picdefs.py

# from enum import Enum
# Code-type
#class code(Enum):
# none = 0
# hipace = 1
# osiris = 2

# Code-type
class code:
  hipace = 'hipace'
  osiris = 'osiris'

# HiPACE definitions   
class hipace:
  # HDF5 definitions  
  class h5:
    # HDF5 GRID dataset
    class g3dkeys:
      beam_charge = 'beam_charge'
    # ...  

    # HDF5 RAW dataset keys
    class rawkeys:
      x1 = 'x1'
      x2 = 'x2'
      x3 = 'x3'
      q = 'q'
      p1 = 'p1'
      p2 = 'p2'
      p3 = 'p3'  
  
    # HDF5 Attribute Keys  
    class attrkeys:
      nx = 'NX'
      xmin = 'XMIN'
      xmax = 'XMAX'
      time = 'TIME'
      dt = 'DT'
      type = 'TYPE'
      name = 'NAME'

# File extensions
class fexts:
    hdf5 = ['.h5','.hdf5']
    ascii = ['.csv','.txt','.ascii']


