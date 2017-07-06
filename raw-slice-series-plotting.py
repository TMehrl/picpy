#!/usr/bin/env python3
# This script may be executed like this:
# nohup ./raw-slice-series-plotting.py <DATA>/ 1> rss.out 2> rss.err &

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
from h5dat import SliceMoms
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
  
  slm = SliceMoms(args[0])

  Xb0 = np.ones(slm.avgx2[0,:].shape)
  for i in range(0,len(slm.zeta_array)):
    if (slm.avgx2[0,i] != 0):
      Xb0[i] = slm.avgx2[0,i]
      
      
  Xb_norm = np.zeros( slm.avgx2.shape )
  zeta_hseed = 1.0
  idx_hseed = (np.abs(slm.zeta_array-zeta_hseed)).argmin()
  
  for i in range(0,len(slm.zeta_array)):
    if (slm.zeta_array[i] <= zeta_hseed): 
      Xb_norm[:,i] = np.absolute( ( slm.avgx2[:,i] - slm.avgx2[:,idx_hseed])/Xb0[i] )

  fig1 = plt.figure()
  cax = plt.pcolormesh( slm.avgx2 )
  cbar = fig1.colorbar(cax)
  cbar.ax.set_ylabel('$X_b$') 
  fig1.savefig(  'Xb_raw.png', 
                format='png')

  fig2 = plt.figure()
  cax = plt.pcolormesh( Xb_norm )
  cbar = fig2.colorbar(cax)
  cbar.ax.set_ylabel('$|X_b/X_{b,0}|$') 
  fig2.savefig(  'Xb.png', 
                format='png')

  fig3 = plt.figure()  
  plt.plot(slm.zeta_array, slm.avgx2[0,:])
  fig3.savefig(  './Xb0.png', 
                format='png')  


  fig4 = plt.figure()
  cax = plt.pcolormesh( np.sqrt( np.absolute( slm.avgx2sq ) ) )
  cbar = fig4.colorbar( cax )
  cbar.ax.set_ylabel('$\sigma_x$') 
  fig4.savefig(  'sigma_x.png', 
                format='png')


  fig5 = plt.figure()
  cax = plt.plot( slm.time_array, np.sqrt( np.absolute( slm.avgx2sq[:,80] )) )
  cbar.ax.set_ylabel('$\sigma_x$') 
  fig5.savefig(  'sigma_x_slice.png', 
                  format='png')

if __name__ == "__main__":
    main()