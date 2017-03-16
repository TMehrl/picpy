#!/usr/bin/env python3

import numpy as np
from optparse import OptionParser
from optparse import OptionGroup
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import picdefs
from grid import Grid3d

def ps_parseopts():

  usg = "Usage: %prog [options] <file or path>"
  
  desc="""This is the picpy postprocessing tool."""
  
  parser = OptionParser(usage=usg, description=desc)
    
  parser.add_option(  "-s", "--save-path", 
                      dest="savepath",
                      metavar="PATH",
                      default='./',
                      help= """Path to which generated files will be saved.
                            (Default: './')""")
  parser.add_option(  "-c", "--code", 
                      type='choice',
                      action='store',
                      dest="piccode",
                      metavar="CODE",
                      choices=['HiPACE', 'OSIRIS',],
                      default='HiPACE',
                      help= """PIC code which was used to generate files
                            (Default: HiPACE).""")
  parser.add_option(  "-d", "--dim", 
                      type='choice',
                      action='store',
                      dest="dimensionality",
                      metavar="DIM",
                      choices=[1, 2, 3,],
                      default=3,
                      help= """Dimensionality of PIC simulation
                            (Default: 3).""")                                     
  parser.add_option(  "-f", "--format", 
                      type='choice',
                      action='store',
                      dest="file_format",
                      metavar="FORMAT",
                      choices=['hdf5', 'bin', 'ascii',],
                      default='hdf5',
                      help= """Format of file (Default: hdf5).""") 
  parser.add_option(  "-a", "--all", 
                      action='store_true',
                      dest="process_all_flag",
                      metavar="DIM",
                      default=False,
                      help="Process all files in path.") 
#   group = OptionGroup(parser, "Options for beam-phase-space (RAW) files",
#                       "These are options for beam-phase-space (RAW) files")
#   group.add_option("-g", action="store_true", help="Group option.")
#   parser.add_option_group(group)

  group = OptionGroup(parser, "Options for grid files",
                      "These are options for grid files")
  group.add_option("-g", action="store_true", help="Group option.")
  parser.add_option_group(group)
  
  return parser



def main():
  
  dir(picdefs)
  
  parser = ps_parseopts()

  (options, args) = parser.parse_args()
  
  if len(args) < 1:
    parser.error("This script requires a file or a path as argument!")
  

  # File 1
  g3d = Grid3d(picdefs.code.HIPACE)
  g3d.read(args[0])
  
  x3_centr_index = math.floor(g3d.nx[2]/2)
  
  fig_zeta = plt.figure()
  plt.pcolormesh( g3d.get_zeta_arr(), 
                  g3d.get_x_arr(1), 
                  g3d.data[:,:,x3_centr_index].T, 
                  cmap='PuBu_r')
  ax = plt.gca()
  ax.set_ylabel(r'$k_p x$', fontsize=14)
  ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
  fig_zeta.savefig('plot_x_zeta.png', format='png')
  
  fig_xi = plt.figure()
  plt.pcolormesh( g3d.get_xi_arr(), 
                  g3d.get_x_arr(1), 
                  g3d.data[:,:,x3_centr_index].T, 
                  cmap='PuBu_r')
  ax = plt.gca()
  ax.set_ylabel(r'$k_p x$', fontsize=14)
  ax.set_xlabel(r'$k_p \xi$', fontsize=14)
  fig_xi.savefig('plot_x_xi.png', format='png')
  
if __name__ == "__main__":
    main()