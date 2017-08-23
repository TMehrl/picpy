#!/usr/bin/env python3
# This script may be executed like this:
# nohup raw-slice-series.py <DATA> 1> rss.out 2> rss.err &

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
from h5dat import HiRAW
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
  indentstr = 'raw_driver_'
  mom_order = 2
  h5fileName = 'slice-avgs.h5'

  dir = DIR(args[0])
  dir.list_files(indentstr)

  if (opts.Nfiles == 0) & (dir.nf > 0):
    Nfiles = dir.nf
  elif int(opts.Nfiles) <= dir.nf:
    Nfiles = int(opts.Nfiles)
  elif dir.nf == 0:
    sys.stderr('Error: No phase space (raw) files in directory!')
    sys.exit() 
  else:
    sys.stderr('Error: Nfiles cannot be smaller than the actual number of files!')

  sys.stdout.write('There are %i raw files to process...\n' % Nfiles)
  sys.stdout.flush()

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
    sys.stdout.write('Processing: %s\t(%i/%i)\n' % (dir.filepath(i), i+1, Nfiles))
    sys.stdout.flush()
    
    raw = HiRAW(dir.filepath(i))
    raw.read_attrs()
    raw.read_data()
    
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
  
  sys.stdout.write('Saving to file: %s\n' % h5fileName)
  sys.stdout.flush()

  h5f = h5py.File(h5fileName, "w")
  dset_zeta_array = h5f.create_dataset( "zeta_array", data = zeta_array )
  dset_time_array = h5f.create_dataset( "time_array", data = time_array )
  dset_avgx1 = h5f.create_dataset(  "avgx1", data = avgx1 )
  dset_avgx2 = h5f.create_dataset(  "avgx2", data = avgx2 )
  dset_avgx3 = h5f.create_dataset(  "avgx3", data = avgx3 )
  dset_avgp1 = h5f.create_dataset(  "avgp1", data = avgp1 )
  dset_avgp2 = h5f.create_dataset(  "avgp2", data = avgp2 )
  dset_avgp3 = h5f.create_dataset(  "avgp3", data = avgp3 )

  if mom_order>1:
    dset_avgx1sq = h5f.create_dataset(  "avgx1sq", data = avgx1sq )
    dset_avgx2sq = h5f.create_dataset(  "avgx2sq", data = avgx2sq )
    dset_avgx3sq = h5f.create_dataset(  "avgx3sq", data = avgx3sq )
    dset_avgp1sq = h5f.create_dataset(  "avgp1sq", data = avgp1sq )
    dset_avgp2sq = h5f.create_dataset(  "avgp2sq", data = avgp2sq )
    dset_avgp3sq = h5f.create_dataset(  "avgp3sq", data = avgp3sq )
    dset_avgx1p1 = h5f.create_dataset(  "avgx1p1", data = avgx1p1 )
    dset_avgx2p2 = h5f.create_dataset(  "avgx2p2", data = avgx2p2 )
    dset_avgx3p3 = h5f.create_dataset(  "avgx3p3", data = avgx3p3 )
      
  h5f.close() 
                                              
  sys.stdout.write('Done!\n')
  sys.stdout.flush()
  
if __name__ == "__main__":
    main()