#!/usr/bin/env python3

import numpy as np
import os
from optparse import OptionParser
from optparse import OptionGroup
import math
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import picdefs
from grid import Grid3d

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
  parser.add_option(  "-c", "--code", 
                      type='choice',
                      action='store',
                      dest="piccode",
                      metavar="CODE",
                      choices = [picdefs.codenames.hipace, picdefs.codenames.osiris,],
                      default = picdefs.codenames.hipace,
                      help= "PIC code which was used to generate files (Default: " +
                            picdefs.codenames.hipace + ").")
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
  parser.add_option(  "-p", "--plane", 
                      type='choice',
                      action='store',
                      dest="plane",
                      metavar="PLANE",
                      choices=['z-x', 'z-y', 'both'],
                      default='z-x',
                      help= """Plane to be plotted (Default: z-x).""")
  parser.add_option(  "-z", "--z-axis", 
                      type='choice',
                      action='store',
                      dest="zax",
                      metavar="ZAXIS",
                      choices=['zeta', 'z', 'xi',],
                      default='zeta',
                      help= """z-axis type (Default: zeta).""")                                                           
  parser.add_option(  "-a", "--all", 
                      action='store_true',
                      dest="process_all",
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


def used_code(code_str):
  if code_str == picdefs.codenames.hipace:
    return picdefs.code.hipace
  elif code_str == picdefs.codenames.osiris:
    return picdefs.code.osiris
    print('Warning: OSIRIS not yet implemented/tested!')
  else:
    print('Error: No/wrong code selected!')
    sys.exit()

def plotfile(file, opts):
  
  code = used_code(opts.piccode)
  
  # File
  g3d = Grid3d(code)
  if opts.verbose == True:  print('Reading: ', file)
  g3d.read(file)
  if opts.verbose == True:  print('Read-in completed.')
  
  if opts.verbose == True: g3d.print_attributes(file)
    
  if opts.zax == 'zeta':
    x_array = g3d.get_zeta_arr()
    xlabel = r'$k_p \zeta$'
  elif opts.zax == 'z':
    x_array = g3d.get_z_arr()
    xlabel = r'$k_p z$'
  elif opts.zax == 'xi':
    x_array = g3d.get_xi_arr()
    xlabel = r'$k_p \xi$'
  else:
    print('Error: No/wrong z-axis option selected!')
    sys.exit()
  
  if opts.plane == 'z-x':
    centr_index = math.floor(g3d.nx[2]/2)
    plane = g3d.data[:,:,centr_index].transpose(1, 0)
    y_array = g3d.get_x_arr(1)
    ylabel = r'$k_p x$'
  elif opts.plane == 'z-y':
    centr_index = math.floor(g3d.nx[1]/2)
    plane = g3d.data[:,centr_index,:].transpose(2, 0, 1)
    y_array = g3d.get_x_arr(2)
    ylabel = r'$k_p y$'
  elif opts.plane == 'both':
    print('Error: Plotting both planes is not implemented yet!')
    sys.exit()
  else:
    print('Error: Wong plane setting!')
    sys.exit()
  
  
  fig = plt.figure()
  plt.pcolormesh( x_array, 
                  y_array, 
                  plane, 
                  cmap='PuBu_r')
  ax = plt.gca()
  ax.set_ylabel(ylabel, fontsize=14)
  ax.set_xlabel(xlabel, fontsize=14)
  fig.savefig(opts.savepath + '/' + 'plot.png', format='png')
  if opts.verbose == True: print('Plot(s) saved at: ' + opts.savepath)    



def main():
  
  parser = ps_parseopts()

  (opts, args) = parser.parse_args()
  
  
  if len(args) < 1:
    parser.error("This script requires a file or a path as argument!")
  
  if os.path.isfile(args[0]):
    file = args[0]
    plotfile(file, opts)
  elif os.path.isdir(args[0]):
    print('Path is dir...!')
    if opts.process_all == True:
      print('Error: Not yet implemented!')
    else:
      print('Error: Use the flag "-a" to process all files in the provided directory!')
      sys.exit() 
  else:
    print('Error: Provided path is neither a file nor a directory!')
    sys.exit()


   
  
if __name__ == "__main__":
    main()