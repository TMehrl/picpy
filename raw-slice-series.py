#!/usr/bin/env python3

import numpy as np
import os
from optparse import OptionParser
from optparse import OptionGroup
import math
import sys
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import picdefs
from h5dat import RAW
from h5dat import DIR
import ps_ana
import h5py

# Parse defaults/definitions
class parsedefs:
  class file_format:
    png = 'png'
    eps = 'eps'
    pdf = 'pdf'
  class save_prefix:
    name = 'g3d_name'


def ps_parseopts():

  usg = "Usage: %prog [options] <file or path>"
  
  desc="""This is the picpy postprocessing tool."""
  
  parser = OptionParser(usage=usg, description=desc)
  parser.add_option(  "-v", "--verbose",
                      action="store_true", 
                      dest="verbose", 
                      default=True,
                      help = "Print info (Default).")
  parser.add_option(  "-q", "--quiet",
                      action="store_false", 
                      dest="verbose",
                      help = "Don't print info.")
  parser.add_option(  "-s", "--save-path", 
                      dest="savepath",
                      metavar="PATH",
                      default='./',
                      help = """Path to which generated files will be saved.
                            (Default: './')""")
  parser.add_option(  "-n", "--name-prefix", 
                      dest="save_prefix",
                      metavar="NAME",
                      default=parsedefs.save_prefix.name,
                      help = """Define customized prefix of output filename.""")  
  parser.add_option(  "-c", "--code", 
                      type='choice',
                      action='store',
                      dest="piccode",
                      metavar="CODE",
                      choices = [picdefs.code.hipace, picdefs.code.osiris,],
                      default = picdefs.code.hipace,
                      help= "PIC code which was used to generate files (Default: " +
                            picdefs.code.hipace + ").")
  parser.add_option(  "-d", "--dim", 
                      type='choice',
                      action='store',
                      dest="dimensionality",
                      metavar="DIM",
                      choices=[1, 2, 3,],
                      default=3,
                      help= """Dimensionality of PIC simulation
                            (Default: 3).""")                                     
  parser.add_option(  "-N", "--number-of-files", 
                      action='store',
                      dest="Nfiles",
                      metavar="NFILES",
                      default=0,
                      help= """Number of files to analyze.""")
#   group = OptionGroup(parser, "Options for beam-phase-space (RAW) files",
#                       "These are options for beam-phase-space (RAW) files")
#   group.add_option("-g", action="store_true", help="Group option.")
#   parser.add_option_group(group)

#   group = OptionGroup(parser, "Options for grid files",
#                       "These are options for grid files")
#   group.add_option("-g", action="store_true", help="Group option.")
#   parser.add_option_group(group)
  
  return parser





def main():
  
  parser = ps_parseopts()

  (opts, args) = parser.parse_args()
  

  nbins = 256
  indentstr = 'raw_beam_phasespace_'
  mom_order = 2
  h5fileName = 'slice-avgs.h5'

  dir = DIR(args[0])
  dir.list_files(indentstr)

  if opts.Nfiles == 0:
    Nfiles = dir.nf
  elif int(opts.Nfiles) <= dir.nf:
    Nfiles = int(opts.Nfiles)
  else:
    print('Error: Nfiles cannot be smaller than the actual number of files!')

  time_array = np.zeros(Nfiles, dtype=np.float32)
  avgx1 = np.zeros((Nfiles, nbins), dtype=np.float32)
  avgx2 = np.zeros((Nfiles, nbins), dtype=np.float32)
  avgx3 = np.zeros((Nfiles, nbins), dtype=np.float32)    
  avgp1 = np.zeros((Nfiles, nbins), dtype=np.float32)
  avgp2 = np.zeros((Nfiles, nbins), dtype=np.float32)
  avgp3 = np.zeros((Nfiles, nbins), dtype=np.float32)
  
  if mom_order>1:
    avgx1sq = np.zeros((Nfiles, nbins), dtype=np.float32)
    avgx2sq = np.zeros((Nfiles, nbins), dtype=np.float32)
    avgx3sq = np.zeros((Nfiles, nbins), dtype=np.float32)    
    avgp1sq = np.zeros((Nfiles, nbins), dtype=np.float32)
    avgp2sq = np.zeros((Nfiles, nbins), dtype=np.float32)
    avgp3sq = np.zeros((Nfiles, nbins), dtype=np.float32)  
    avgx1p1 = np.zeros((Nfiles, nbins), dtype=np.float32)
    avgx2p2 = np.zeros((Nfiles, nbins), dtype=np.float32)
    avgx3p3 = np.zeros((Nfiles, nbins), dtype=np.float32)
         
  for i in range(0, Nfiles):
    print('Processing: %s' % dir.filepath(i))
    raw = RAW(dir.filepath(i), picdefs.code.hipace)
    time_array[i] = raw.time
    slices = ps_ana.SLICES(raw, nbins=nbins)
    slices.calc_moments(order = mom_order)
    avgx1[i,:] = slices.avgx1      
    avgx2[i,:] = slices.avgx2
    avgx3[i,:] = slices.avgx3
    avgp1[i,:] = slices.avgp1      
    avgp2[i,:] = slices.avgp2
    avgp3[i,:] = slices.avgp3  
    
    if mom_order>1:
      avgx1sq[i,:] = slices.avgx1sq      
      avgx2sq[i,:] = slices.avgx2sq
      avgx3sq[i,:] = slices.avgx3sq
      avgp1sq[i,:] = slices.avgp1sq    
      avgp2sq[i,:] = slices.avgp2sq
      avgp3sq[i,:] = slices.avgp3sq
      avgx1p1[i,:] = slices.avgx1p1    
      avgx2p2[i,:] = slices.avgx2p2
      avgx3p3[i,:] = slices.avgx3p3

  zeta_array = slices.centers

  h5f = h5py.File(h5fileName, "w")
  dset = h5f.create_dataset("zeta_array", zeta_array.shape, zeta_array.dtype)
  dset = h5f.create_dataset("time_array", time_array.shape, time_array.dtype)
  dset = h5f.create_dataset("avgx1", avgx1.shape, avgx1.dtype)
  dset = h5f.create_dataset("avgx2", avgx2.shape, avgx2.dtype)
  dset = h5f.create_dataset("avgx3", avgx3.shape, avgx3.dtype)
  dset = h5f.create_dataset("avgp1", avgp1.shape, avgp1.dtype)
  dset = h5f.create_dataset("avgp2", avgp2.shape, avgp2.dtype)
  dset = h5f.create_dataset("avgp3", avgp3.shape, avgp3.dtype)
  if mom_order>1:
    dset = h5f.create_dataset("avgx1sq", avgx1sq.shape, avgx1sq.dtype)
    dset = h5f.create_dataset("avgx2sq", avgx2sq.shape, avgx2sq.dtype)
    dset = h5f.create_dataset("avgx3sq", avgx3sq.shape, avgx3sq.dtype)
    dset = h5f.create_dataset("avgp1sq", avgp1sq.shape, avgp1sq.dtype)
    dset = h5f.create_dataset("avgp2sq", avgp2sq.shape, avgp2sq.dtype)
    dset = h5f.create_dataset("avgp3sq", avgp3sq.shape, avgp3sq.dtype)   
    dset = h5f.create_dataset("avgx1p1", avgx1p1.shape, avgx1p1.dtype)
    dset = h5f.create_dataset("avgx2p2", avgx2p2.shape, avgx2p2.dtype)
    dset = h5f.create_dataset("avgx3p3", avgx3p3.shape, avgx3p3.dtype) 
  h5f.close() 
               
  fig = plt.figure()  
  plt.plot(slices.centers, slices.avgx2)
  fig.savefig(  './Xb.png', 
                 format='png') 

  fig = plt.figure()
  cax = plt.pcolormesh( avgx2 )
  fig.savefig(  './Xb_evolv.png', 
                 format='png')                               
  
if __name__ == "__main__":
    main()