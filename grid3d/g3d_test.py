#!/usr/bin/env python3

import h5py
import numpy as np
import sys
from scipy import constants
import os
import csv
import math
from enum import Enum
from optparse import OptionParser
import matplotlib.pyplot as plt

# Global
h5_ext_list = ['.h5','.hdf5']
ascii_ext_list = ['.csv','.txt','.ascii']
xlim_attr_str = ['XMIN','XMAX']

class PIC_CODE(Enum):
  NONE = 0
  HIPACE = 1
  OSIRIS = 2


class GRID3D:
  def __init__(self, file):
    fpath, fext = os.path.splitext(file)
    self.filename = file
    with h5py.File(file,'r') as hf:
      self.attrs = hf.attrs
      self.h5keys = list(hf.keys())
      if len(self.h5keys) == 0:
        print('ERROR:\tHDF5 file "%s" does not contain any dataset!' %file)
        sys.exit()
      elif len(self.h5keys) == 1:
        self.data = np.array(hf.get(self.h5keys[0]))
      else:
        print('ERROR:\tHDF5 file "%s" contains more than one dataset!' %file)
        sys.exit()
      for item in self.attrs.keys():
        print(item + ":", self.attrs[item])
      hf.close()
#  def get_x_arr(self,dim):

def read_file(fname):

  fpath, fext = os.path.splitext(fname)

  print('Reading', fname)
  if any(fext == s for s in h5_ext_list):
    grid3d = GRID3D(fname)
  else:
    print('ERROR:\tExtension of file "%s" is not supported!' %fname)
    print('\tAllowed hdf5-file extensions: ',list(h5_ext_list))
    sys.exit()
  print('Reading complete.')
  return grid3d
  
  
  
  
def main():

  NUM_ARGS = 1

  usage = "usage: %prog [options] <file>"
  parser = OptionParser(usage=usage)
  (options, args) = parser.parse_args()

  if len(args) != NUM_ARGS:
    parser.error("This script requires an argument!")  

  # File 1
  g3d = read_file(args[0])
  
#  for x in range(0, 3):
#    print('We are on time %d' % (x))
  
  
#  if xlim_attr_str[0] in g3d.attrs.keys():
#    print('Found: %s',xlim_attr_str[0])

  plt.figure()
  plt.pcolormesh(g3d.data[:,:,64])
  plt.show()
  
if __name__ == "__main__":
    main()