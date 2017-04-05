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
                      metavar="FORMAT",
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
  

  dir = DIR(args[0])
  dir.list_files('raw_beam_phasespace_')

  nbins=256

  if opts.Nfiles == 0:
    Nfiles = dir.nf
  elif int(opts.Nfiles) <= dir.nf:
    Nfiles = int(opts.Nfiles)
  else:
    print('Error: Nfiles cannot be smaller than the actual number of files!')

  t_array = np.zeros(Nfiles, dtype=np.float32)
  AVGx2 = np.zeros((Nfiles, nbins), dtype=np.float32)
  AVGx3 = np.zeros((Nfiles, nbins), dtype=np.float32)



  for i in range(0,Nfiles):
    raw = RAW(dir.filepath(i), picdefs.code.hipace)
    t_array[i] = raw.time
    print(dir.filepath(i))
    slices = ps_ana.SLICES(raw, nbins=nbins)
    slices.calc_moments()
      
    AVGx2[i,:] = slices.avgx2
    AVGx3[i,:] = slices.avgx3
  
    
  fig = plt.figure()  
  plt.plot(slices.centers, slices.avgx2)
  fig.savefig(  './Xb.png', 
                 format='png') 

  fig = plt.figure()
  cax = plt.pcolormesh( AVGx2 )
  fig.savefig(  './Xb_evolv.png', 
                 format='png') 
#   # File
#   raw = RAW(file, picdefs.code.hipace)
#   raw.print_datasets(file)
#   raw.print_attributes(file)
#   
#   
#   
#   fig = plt.figure()
#   plt.hist2d(raw.x1, raw.p1, bins=[256, 512], norm=LogNorm())
#   plt.colorbar()
#   ax = plt.gca()
#   ax.set_ylabel('p1', fontsize=14)
#   ax.set_xlabel('x1', fontsize=14)
#   fig.savefig(  './long_ps.png', 
#                 format='png')
#   
#   slices = ps_ana.SLICES(raw, nbins=256)
#   slices.calc_moments()
#   fig = plt.figure()  
#   plt.plot(slices.centers, np.sqrt(slices.avgx2sq))
#   fig.savefig(  './sigma_x.png', 
#                 format='png')
#                 
#   fig = plt.figure()  
#   plt.plot(slices.centers, slices.charge)
#   fig.savefig(  './curr_profile.png', 
#                 format='png') 
#                 
#   fig = plt.figure()  
#   plt.plot(slices.centers, slices.avgx2)
#   fig.savefig(  './Xb.png', 
#                 format='png')                                
  
if __name__ == "__main__":
    main()