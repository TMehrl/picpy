#!/usr/bin/env python3
# picdefs.py

from enum import Enum

# Code-type
class code(Enum):
 none = 0
 hipace = 1
 osiris = 2

class codenames:
  hipace = 'hipace'
  osiris = 'osiris'

# HiPACE definitions   
class hipace:
  # HDF5 Grid 3D Keys  
  class h5attrkeys:
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


