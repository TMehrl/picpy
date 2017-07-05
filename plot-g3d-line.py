#!/usr/bin/env python3

import numpy as np
import os
import argparse
import math
import sys
import matplotlib
# Force matplotlib to not use any Xwindows backend.
# matplotlib.use('Agg')
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
    xy = 'xy'    
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


def two_floats(value):
    values = value.split()
    if len(values) != 2:
        raise argparse.ArgumentError
    values = map(float, values)
    return values


def ps_parseargs():

  usg = "Usage: %prog [options] <file or path>"
  
  desc="""This is the picpy postprocessing tool."""
  
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument(  'path',
                        nargs="+",
                        metavar='PATH', 
                        help='Path to file.')
  parser.add_argument(  "-v", "--verbose",
                        dest="verbose", 
                        action="store_true", 
                        default=True,
                        help = "Print info (Default).")
  parser.add_argument(  "--show",
                        dest="ifshow", 
                        action="store_true", 
                        default=False,
                        help = "Show figure.")
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
                        choices=[ parsedefs.plane.zx, parsedefs.plane.zy,
                                  parsedefs.plane.xy,],
                        default=parsedefs.plane.zx,
                        help= """Plane to be plotted (Default: zx).""")  
  parser.add_argument(  "--plane-index", 
                        action='store',
                        dest="plane_index",
                        metavar="PLANE-IDX",
                        default=None,
                        type=int,
                        help= """Index of plane.""")                                              
  parser.add_argument(  "-z", "--z-axis", 
                        action='store',
                        dest="zax",
                        metavar="ZAXIS",
                        choices=[ parsedefs.zax.zeta, 
                                parsedefs.zax.z, 
                                parsedefs.zax.xi,],
                        default=parsedefs.zax.zeta,
                        help= "z-axis type (Default: " + parsedefs.zax.zeta + ").")
  parser.add_argument(  "--cscale",
                        action='store',
                        dest="cscale",
                        metavar="CSCALE",
                        choices=[ "lin", "log",],
                        default="lin",
                        help= "z-axis type (Default: " + parsedefs.zax.zeta + ").")                                                                                  
  parser.add_argument(  "-a", "--all", 
                        action='store_true',
                        dest="process_all",
                        default=False,
                        help="Process all files in path.")  
  parser.add_argument(  '--cblim', 
                        help='Colorbar axis limits',
                        action='store', 
                        dest="cblim",
                        metavar="'CBMIN CBMAX'",                        
                        type=two_floats,
                        default=None)
                     
  
  return parser



def plotfile( args ):


  fig = plt.figure()
  savename = ""
  
  for file in args.path: 

    # File
    if args.verbose == True:  print('Reading: ', file)
    g3d = Grid3d(file, args.piccode)
    if args.verbose == True:  print('Read-in completed.')
  
    if args.verbose == True: 
      g3d.print_datasets()
      g3d.print_attributes()
    
    if (args.plane == parsedefs.plane.zx) | (args.plane == parsedefs.plane.zy):
    
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

    x1_index = math.floor(g3d.nx[1]/2) - 1
    x2_index = math.floor(g3d.nx[2]/2) - 1
    
    data = g3d.data[:,x1_index,x2_index]

    saveformat = args.file_format  
    filesuffix = '_%06.f' % (np.floor(g3d.time))
  
    if args.save_prefix != parsedefs.save_prefix.name:
      fileprefix = args.save_prefix
    else: 
      fileprefix = g3d.name
  
    print(type(args.cblim))  
  
    savename += fileprefix + filesuffix + '.'

    if args.cscale == "log":
      data = np.log(abs(data))
    
    plt.plot( x_array, 
              data,
              label=file)
    ax = plt.gca()
    ax.set_ylabel(g3d.name, fontsize=14)
    ax.set_xlabel(xlabel, fontsize=14)
  
  
  plt.legend()
  savename += saveformat
    
  fig.savefig(  args.savepath + '/' + savename, 
                format=saveformat)
  if args.verbose: print('Saved "' + savename + '" at: ' + args.savepath)    

  if args.ifshow: plt.show()



def main():
  
  parser = ps_parseargs()

  args = parser.parse_args()
    
  if all( os.path.isfile(str) for str in args.path ):
    plotfile( args )
  
  elif ( len(args.path)==1 ) & ( os.path.isdir( args.path ) ):
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