#!/usr/bin/env python3

import numpy as np
import os
from optparse import OptionParser
from optparse import OptionGroup
import math
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import picdefs
from h5dat import RAW
import ps_ana


def main():
  
  usg = "Usage: %prog [options] <file or path>"
  
  desc="""This is the picpy postprocessing tool."""
  
  parser = OptionParser(usage=usg, description=desc)
  (opts, args) = parser.parse_args()
    
  if len(args) < 1:
    parser.error("This script requires a file or a path as argument!")
  
  if os.path.isfile(args[0]):
    file = args[0]
  
  # File
  raw = RAW(file, picdefs.code.hipace)
  raw.print_datasets(file)
  raw.print_attributes(file)
  
  # Force matplotlib to not use any Xwindows backend.
  matplotlib.use('Agg')
  
  fig = plt.figure()
  plt.hist2d(raw.x1, raw.p1, bins=[256, 512], norm=LogNorm())
  plt.colorbar()
  ax = plt.gca()
  ax.set_ylabel('p1', fontsize=14)
  ax.set_xlabel('x1', fontsize=14)
  fig.savefig(  './long_ps.png', 
                format='png')


  moments = ps_ana.moments(raw.x2, raw.p2, raw.q, order=2, central=False)
  print(moments)
  
  slices = ps_ana.SLICES(raw)
  slices.calc_moments()
  fig = plt.figure()  
  plt.plot(slices.centers, np.sqrt(slices.x2p2_moms[:,0]))
  fig.savefig(  './sigma_x.png', 
                format='png')
  
if __name__ == "__main__":
    main()