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
  parser.add_option(  "-n", "--name-prefix", 
                      dest="save_name_prefix",
                      metavar="NAME",
                      default='g3d_slice',
                      help = """Prefix of output filename.
                            (Default: 'g3d_slice')""")  
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
  parser.add_option(  "-f", "--format", 
                      type='choice',
                      action='store',
                      dest="file_format",
                      metavar="FORMAT",
                      choices=['png', 'pdf', 'eps',],
                      default='png',
                      help= """Format of output file (Default: png).""")
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

#   group = OptionGroup(parser, "Options for grid files",
#                       "These are options for grid files")
#   group.add_option("-g", action="store_true", help="Group option.")
#   parser.add_option_group(group)
  
  return parser



def plotfile(file, opts):
    
  # File
  g3d = Grid3d(opts.piccode)
  
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
    data = g3d.data[:,:,centr_index].transpose(1, 0)
    y_array = g3d.get_x_arr(1)
    ylabel = r'$k_p x$'
  elif opts.plane == 'z-y':
    centr_index = math.floor(g3d.nx[1]/2)
    data = g3d.data[:,centr_index,:].transpose(2, 0, 1)
    y_array = g3d.get_x_arr(2)
    ylabel = r'$k_p y$'
  elif opts.plane == 'both':
    print('Error: Plotting both planes is not implemented yet!')
    sys.exit()
  else:
    print('Error: Wong plane setting!')
    sys.exit()

  saveformat = opts.file_format  
  filesuffix = '_%06.f' % (g3d.time)
  savename = opts.save_name_prefix + filesuffix + '.' + saveformat 

  
  fig = plt.figure()
  plt.pcolormesh( x_array, 
                  y_array, 
                  data, 
                  cmap='PuBu_r')
  ax = plt.gca()
  ax.set_ylabel(ylabel, fontsize=14)
  ax.set_xlabel(xlabel, fontsize=14)
  fig.savefig(  opts.savepath + '/' + savename, 
                format=saveformat)
  if opts.verbose == True: print('Saved "' + savename + '" at: ' + opts.savepath)    


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
      sys.exit()
    else:
      print('Error: Use the flag "-a" to process all files in the provided directory!')
      sys.exit() 
  
  else:
    print('Error: Provided path is neither a file nor a directory!')
    sys.exit()


   
  
if __name__ == "__main__":
    main()