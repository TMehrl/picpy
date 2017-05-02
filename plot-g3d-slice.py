#!/usr/bin/env python3

import numpy as np
import os
from argparse import ArgumentParser
import math
import sys
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
import picdefs
from h5dat import Grid3d


# Parse defaults/definitions
class parsedefs:
  class plane:
    zx = 'zx'
    zy = 'zy'
  class file_format:
    png = 'png'
    eps = 'eps'
    pdf = 'pdf'
  class zax:
    zeta  = 'zeta'
    z     = 'z'
    xi    = 'xi'
  class save_prefix:
    name = 'g3d_name'



def ps_parseargs():

  usg = "Usage: %prog [options] <file or path>"
  
  desc="""This is the picpy postprocessing tool."""
  
  parser = ArgumentParser(description=desc)
  parser.add_argument(  'path', 
                        metavar='PATH', 
                        help='Path to file.')
  parser.add_argument('--sum', dest='accumulate', action='store_const',
                      const=sum, default=max,
                      help='sum the integers (default: find the max)')
  parser.add_argument(  "-v", "--verbose",
                        dest="verbose", 
                        action="store_true", 
                        default=True,
                        help = "Print info (Default).")
  parser.add_argument(  "-q", "--quiet",
                        dest="verbose",
                        action="store_false", 
                        help = "Don't print info.")
  parser.add_argument(  "-s", "--save-path", 
                        action="store", 
                        dest="savepath",
                        metavar="PATH",
                        default='./',
                        help = """Path to which generated files will be saved.
                            (Default: './')""")
  parser.add_argument(  "-n", "--name-prefix", 
                        action="store", 
                        dest="save_prefix",
                        metavar="NAME",
                        default=parsedefs.save_prefix.name,
                        help = """Define customized prefix of output filename.""")  
  parser.add_argument(  "-c", "--code", 
                        action="store", 
                        dest="piccode",
                        metavar="CODE",
                        choices = [picdefs.code.hipace, picdefs.code.osiris,],
                        default = picdefs.code.hipace,
                        help="PIC code (Default: " + picdefs.code.hipace + ").")
  parser.add_argument(  "-d", "--dim", 
                        action='store',
                        dest="dimensionality",
                        metavar="DIM",
                        choices=[1, 2, 3,],
                        default=3,
                        help= """Dimensionality of PIC simulation
                            (Default: 3).""")                                     
  parser.add_argument(  "-f", "--format", 
                        action='store',
                        dest="file_format",
                        metavar="FORMAT",
                        choices=[ parsedefs.file_format.png, 
                                  parsedefs.file_format.pdf, 
                                  parsedefs.file_format.eps,],
                        default=parsedefs.file_format.png,
                        help= """Format of output file (Default: png).""")
  parser.add_argument(  "-p", "--plane", 
                        action='store',
                        dest="plane",
                        metavar="PLANE",
                        choices=[ parsedefs.plane.zx, parsedefs.plane.zy,],
                        default=parsedefs.plane.zx,
                        help= """Plane to be plotted (Default: zx).""")
  parser.add_argument(  "-z", "--z-axis", 
                        action='store',
                        dest="zax",
                        metavar="ZAXIS",
                        choices=[ parsedefs.zax.zeta, 
                                parsedefs.zax.z, 
                                parsedefs.zax.xi,],
                        default=parsedefs.zax.zeta,
                        help= "z-axis type (Default: " + parsedefs.zax.zeta + ").")                                                           
  parser.add_argument(  "-a", "--all", 
                        action='store_true',
                        dest="process_all",
                        default=False,
                        help="Process all files in path.")

#   parser.add_argument(  '--cblim', 
#                       help = 'Colorbar axis limits',
#                       action = 'store', 
#                       type=two_floats,
#                       default=[-1.0, 0.0])
                     
#   group = OptionGroup(parser, "Options for beam-phase-space (RAW) files",
#                       "These are options for beam-phase-space (RAW) files")
#   group.add_argument("-g", action="store_true", help="Group option.")
#   parser.add_argument_group(group)

#   group = OptionGroup(parser, "Options for grid files",
#                       "These are options for grid files")
#   group.add_argument("-g", action="store_true", help="Group option.")
#   parser.add_argument_group(group)
  
  return parser



def plotfile(file, args):
    
  # File
  if args.verbose == True:  print('Reading: ', file)
  g3d = Grid3d(file, args.piccode)
  if args.verbose == True:  print('Read-in completed.')
  
  if args.verbose == True: 
    g3d.print_datasets(file)
    g3d.print_attributes(file)
    
  if args.zax == parsedefs.zax.zeta:
    x_array = g3d.get_zeta_arr()
    xlabel = r'$k_p \zeta$'
  elif args.zax == parsedefs.zax.z:
    x_array = g3d.get_z_arr()
    xlabel = r'$k_p z$'
  elif args.zax == parsedefs.zax.xi:
    x_array = g3d.get_xi_arr()
    xlabel = r'$k_p \xi$'
  else:
    print('Error: No/wrong z-axis option selected!')
    sys.exit()
  
  if args.plane == parsedefs.plane.zx:
    centr_index = math.floor(g3d.nx[2]/2) - 1
    if g3d.nx[2]%2 == 1:
      data = g3d.data[:,:,centr_index].transpose(1, 0)
    else:
      data = ( g3d.data[:,:,centr_index].transpose(1, 0)
                + g3d.data[:,:,centr_index+1].transpose(1, 0) )/2
    y_array = g3d.get_x_arr(1)
    ylabel = r'$k_p x$'
  elif args.plane == parsedefs.plane.zy:
    centr_index = math.floor(g3d.nx[1]/2) - 1
    if g3d.nx[1]%2 == 1:
      data = g3d.data[:,centr_index,:].transpose(1, 0)
    else:
      data = ( g3d.data[:,centr_index,:].transpose(1, 0)
                + g3d.data[:,centr_index+1,:].transpose(1, 0) )/2   
    y_array = g3d.get_x_arr(2)
    ylabel = r'$k_p y$'
  else:
    print('Error: Wong plane setting!')
    sys.exit()

  saveformat = args.file_format  
  filesuffix = '_%06.f' % (np.floor(g3d.time))
  
  if args.save_prefix != parsedefs.save_prefix.name:
    fileprefix = args.save_prefix
  else: 
    fileprefix = g3d.name
  
  
  
  savename = fileprefix + filesuffix + '.' + saveformat
  
  if g3d.type == picdefs.hipace.h5.g3dtypes.density:
    colormap = 'PuBu_r';
    cbmax = np.amax(data)
    cbmin = np.amin(data)
  elif g3d.type == picdefs.hipace.h5.g3dtypes.field:
    colormap = cm.coolwarm
    cbmax = np.amax(abs(data))
    cbmin = -cbmax

  cbmin = -0.001
  cbmax = 0.001

#   if args.cbmin != None:
#     cbmin = args.cbmin
# 
#   if args.cbmax != None:
#     cbmax = args.cbmax

  fig = plt.figure()
  cax = plt.pcolormesh( x_array, 
                        y_array, 
                        data,
                        vmin=cbmin, vmax=cbmax,
                        cmap=colormap)
  ax = plt.gca()
  ax.set_ylabel(ylabel, fontsize=14)
  ax.set_xlabel(xlabel, fontsize=14)
  cbar = fig.colorbar(cax)
  cbar.ax.set_ylabel(g3d.name)
  
  fig.savefig(  args.savepath + '/' + savename, 
                format=saveformat)
  if args.verbose: print('Saved "' + savename + '" at: ' + args.savepath)    


def main():
  
  if len(sys.argv) < 1:
    parser.error("This script requires a file or a path as argument!")  
  
  parser = ps_parseargs()

  args = parser.parse_args()
    
  if os.path.isfile(args.path):
    file = args.path
    plotfile(file, args)
  
  elif os.path.isdir(args.path):
    print('Path is dir...!')
    if args.process_all == True:
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