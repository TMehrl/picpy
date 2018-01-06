#!/usr/bin/env python3

import numpy as np
import os
import argparse
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

def two_ints(value):
    values = value.split()
    if len(values) != 2:
        raise argparse.ArgumentError
    values = map(int, values)
    return values

def parser(ptype='none'):

  usg = "Usage: %prog [options] <file or path>"
  
  desc="""This is the picpy postprocessing tool."""
  
  savepath = './plots'
  file_format = None

  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument(  'path', 
                        metavar='PATH',
                        nargs='*',
                        help='Path to grid file.')
  parser.add_argument(  "-v", "--verbose",
                        dest="verbose", 
                        action="store_true", 
                        default=True,
                        help = "Print info (Default).")
  parser.add_argument(  "-q", "--quiet",
                        dest="verbose",
                        action="store_false", 
                        help = "Don't print info.")
  parser.add_argument(  "--show",
                        dest="ifshow", 
                        action="store_true", 
                        default=False,
                        help = "Show figure.")
  parser.add_argument(  "-a", "--all", 
                        action='store_true',
                        dest="process_all",
                        default=False,
                        help="Process all files in path.") 
  parser.add_argument(  "--name-prefix", 
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
  parser.add_argument(  "-z", "--z-axis", 
                        action='store',
                        dest="zax",
                        metavar="ZAXIS",
                        choices=[ parsedefs.zax.zeta, 
                                parsedefs.zax.z, 
                                parsedefs.zax.xi,],
                        default=parsedefs.zax.zeta,
                        help= "z-axis type (Default: " + parsedefs.zax.zeta + ").")
  if ptype == 'slice':
    # Slice plot specific arguments
    savepath += '/g3d-slice'
    file_format = parsedefs.file_format.png
    parser.add_argument(  "-p", "--plane", 
                          action='store',
                          dest="plane",
                          metavar="PLANE",
                          choices=[ 'xy', 'yz', 'xz',
                                    'yx', 'zy', 'zx'],
                          default='zx',
                          help= """Plane to be plotted (Default: zx).""")  
    parser.add_argument(  "--plane-index", 
                          action='store',
                          dest="plane_index",
                          metavar="PLANE-IDX",
                          default=None,
                          type=int,
                          help= """Index of plane.""")  
    parser.add_argument(  "--cscale",
                          action='store',
                          dest="cscale",
                          metavar="CSCALE",
                          choices=[ "lin", "log",],
                          default="lin",
                          help= "z-axis type (Default: " + parsedefs.zax.zeta + ").")                                                                                   
    parser.add_argument(  '--cblim', 
                          help='Colorbar axis limits',
                          action='store', 
                          dest="cblim",
                          metavar="'CBMIN CBMAX'",                        
                          type=two_floats,
                          default=None)
  elif ptype == 'line':
    # Line plot specific arguments
    savepath += '/g3d-line'
    file_format = parsedefs.file_format.eps
    parser.add_argument(  "-l", "--lineout-axis", 
                          action='store',
                          dest="loutax",
                          metavar="LOUTAX",
                          choices=[ 'x', 'y', 'z'],
                          default='z',
                          help= """Axis along which lineout is generated (Default: z).""")     
    parser.add_argument(  '--lineout-indices', 
                          help='Indices for which lineout is taken.',
                          action='store', 
                          dest="lout_idx",
                          metavar="'idx0 idx1'",                        
                          type=two_ints,
                          default=None)

  # General arguments with plot type-specific defaults
  parser.add_argument(  "-s", "--save-path", 
                        action="store", 
                        dest="savepath",
                        metavar="PATH",
                        default=savepath,
                        help = """Path to which generated files will be saved.
                            (Default: './')""")
  parser.add_argument(  "-f", "--format", 
                        action='store',
                        dest="file_format",
                        metavar="FORMAT",
                        choices=[ parsedefs.file_format.png, 
                                  parsedefs.file_format.pdf, 
                                  parsedefs.file_format.eps,],
                        default=file_format,
                        help= """Format of output file (Default: png).""")                               
  return parser


def gen_pretty_grid_name( gname ):
  if gname == 'ExmBy':
    return r'$E_x-B_y$'
  elif gname == 'EypBx':
    return r'$E_y+B_x$'
  elif gname == 'Ez':
    return r'$E_z$'
  elif gname == 'Bx':
    return r'$B_x$'
  elif gname == 'By':
    return r'$B_y$'
  elif gname == 'Bz':
    return r'$B_z$'
  elif gname == 'Jx':
    return r'$J_x$'
  elif gname == 'Jy':
    return r'$J_y$'
  elif gname == 'Jz':
    return r'$J_z$'    
  elif gname == 'Jz':
    return r'$J_z$'
  elif gname == 'plasma_charge':
    return r'$\rho_p$'
  elif gname == 'beam_charge':
    return r'$\rho_b$'
  else:
    return gname  

def is_h5_file(fext):
  return any(fext == h5ext for h5ext in picdefs.fexts.hdf5) 

def is_g3d_file(fname):
  return any((mq in fname) for mq in picdefs.hipace.h5.g3dtypes.list)

def is_h5g3d_file(file):
  fname, fext = os.path.splitext(file)
  return is_h5_file(fext) and is_g3d_file(fname)

# General Grid3D_plot class
class G3d_plot:
  def __init__(self, file, args):
    self.args = args
    self.file = file  

    # Reading hdf5 attributes:
    if self.args.verbose:  print('Getting attributes of ', file)
    self.g3d = Grid3d(file)
    if self.args.verbose:     
      self.g3d.print_datasets()
      self.g3d.print_attributes()

  def set_xaxis(self, xax_str):
    # define axis labels and arrays
    if xax_str == 'z':
      if self.args.zax == parsedefs.zax.zeta:
        self.x_array = self.g3d.get_zeta_arr()
        self.xlabel = r'$k_p \zeta$'
      elif self.args.zax == parsedefs.zax.z:
        self.x_array = self.g3d.get_z_arr()
        self.xlabel = r'$k_p z$'
      elif self.args.zax == parsedefs.zax.xi:
        self.x_array = self.g3d.get_xi_arr()
        self.xlabel = r'$k_p \xi$'
      else:
        print('Error: No/wrong z-axis option selected!')
        sys.exit()
    elif xax_str == 'x':    
      self.x_array = self.g3d.get_x_arr(1)
      self.xlabel = r'$k_p x$'
    elif xax_str == 'y':
      self.x_array = self.g3d.get_x_arr(2)
      self.xlabel = r'$k_p y$'        
    else:
      print('Error: Wrong x-axis string!')
      sys.exit()

  def mkdirs_if_nexist( self ):
    folders = []
    path = self.args.savepath

    while 1:
      path, folder = os.path.split(path)
      if folder != "":
        folders.append(folder)
      else:
        if path != "":
            folders.append(path)
        break
    folders.reverse()

    path = ""
    for folder in folders:
      path = path + folder + "/"
      if not os.path.isdir(path):
        print("Creating folder: " + path)
        os.mkdir(path)

class G3d_plot_slice(G3d_plot):
  def __init__(self, file, args):
    G3d_plot.__init__(self, file, args)

    self.set_xaxis( self.args.plane[0] )
    self.set_yaxis( self.args.plane[1] )
    if self.args.verbose: print('Reading data...')
    self.read() 
    if self.args.verbose: print('Read-in completed.')
    self.set_cmap()

  def set_yaxis(self, yax_str):
    # define axis labels and arrays
    if yax_str == 'z':
      if self.args.zax == parsedefs.zax.zeta:
        self.y_array = self.g3d.get_zeta_arr()
        self.ylabel = r'$k_p \zeta$'
      elif self.args.zax == parsedefs.zax.z:
        self.y_array = self.g3d.get_z_arr()
        self.ylabel = r'$k_p z$'
      elif self.args.zax == parsedefs.zax.xi:
        self.y_array = self.g3d.get_xi_arr()
        self.ylabel = r'$k_p \xi$'
      else:
        print('Error: No/wrong z-axis option selected!')
        sys.exit()
    elif yax_str == 'x':    
      self.y_array = self.g3d.get_x_arr(1)
      self.ylabel = r'$k_p x$'
    elif yax_str == 'y':
      self.y_array = self.g3d.get_x_arr(2)
      self.ylabel = r'$k_p y$'        
    else:
      print('Error: Wrong y-axis string!')
      sys.exit()

  def read( self ):
    # read slice
    if 'z' in self.args.plane:      
      if 'x' in self.args.plane:
        if self.args.plane_index == None:
          index = math.floor(self.g3d.nx[2]/2) - 1
          self.slice = self.g3d.read_slice(i2=index)
          if self.g3d.nx[2]%2 == 0:
            self.slice = (self.slice + self.g3d.read_slice(i2=index+1))/2
        else:
          self.slice = self.g3d.read_slice(i2=self.args.plane_index)
      elif 'y' in self.args.plane:
        if self.args.plane_index == None:
          index = math.floor(self.g3d.nx[1]/2) - 1
          self.slice = self.g3d.read_slice(i1=index)
          if self.g3d.nx[1]%2 == 0:
            self.slice = ( self.slice  + self.g3d.read_slice(i1=index+1) )/2
        else:
          self.slice = self.g3d.read_slice(i1=self.args.plane_index)
    elif ('x' in self.args.plane) and ('y' in self.args.plane):
      if self.args.plane_index == None:    
        index = math.floor(self.g3d.nx[0]/2) - 1
        self.slice = self.g3d.read_slice(i0=index)   
        if self.g3d.nx[0]%2 == 0:
          self.slice = ( self.slice  + self.g3d.read_slice(i0=index+1) )/2 
      else:
        self.slice = self.g3d.read_slice(i0=self.args.plane_index)

    if self.args.plane in ['xy','zx','zy']:
      self.slice = np.transpose( self.slice )

  def set_cmap( self ):
    if self.args.cscale == "log":
      self.slice = np.log(abs(self.slice))
      
    cblim = [0.0, 0.0]
    
    if self.g3d.type == picdefs.hipace.h5.g3dtypes.density:
      self.colormap = 'PuBu_r';
      cblim[0] = np.amin(self.slice)
      cblim[1] = np.amax(self.slice)
    elif self.g3d.type == picdefs.hipace.h5.g3dtypes.field:
      self.colormap = cm.coolwarm
      cblim[0] = -np.amax(abs(self.slice))
      cblim[1] = np.amax(abs(self.slice))
    elif self.g3d.type == picdefs.hipace.h5.g3dtypes.current:
      self.colormap = cm.coolwarm
      cblim[0] = -np.amax(abs(self.slice))
      cblim[1] = np.amax(abs(self.slice))

    if self.args.cblim != None:
      cblim = list(self.args.cblim)
    
    self.cblim = cblim

  def plot( self, ifsave=True ):  
    if self.args.verbose: print('Generating slice plot') 
    saveformat = self.args.file_format  
    filesuffix = '_%06.f' % (np.floor(self.g3d.time))
    
    if self.args.save_prefix != parsedefs.save_prefix.name:
      fileprefix = self.args.save_prefix
    else: 
      fileprefix = self.g3d.name

    savename = fileprefix + filesuffix + '.' + saveformat

    fig = plt.figure()
    cax = plt.pcolormesh( self.x_array, 
                          self.y_array, 
                          self.slice,
                          vmin=self.cblim[0], vmax=self.cblim[1],
                          cmap=self.colormap)
    ax = plt.gca()
    ax.set_ylabel(self.ylabel, fontsize=14)
    ax.set_xlabel(self.xlabel, fontsize=14)
    cbar = fig.colorbar(cax)
    cbar.ax.set_ylabel( gen_pretty_grid_name( self.g3d.name ), fontsize=14 )

    self.mkdirs_if_nexist()

    if saveformat==parsedefs.file_format.png:
      fig.savefig(  self.args.savepath + '/' + savename, 
                format=saveformat,
                dpi=600)
    else:    
      fig.savefig(  self.args.savepath + '/' + savename, 
                    format=saveformat)
    if self.args.verbose: print('Saved "' + savename + '" at: ' + self.args.savepath)    
    
    if self.args.ifshow: plt.show()
    plt.close(fig)



class G3d_plot_line(G3d_plot):
  def __init__(self, file, args):
    G3d_plot.__init__(self, file, args)

    self.set_xaxis( self.args.loutax )
    if self.args.verbose: print('Reading data...')
    self.read() 
    if self.args.verbose: print('Read-in completed.')
    self.set_yaxis()

  def set_yaxis( self ):
    self.ylabel = gen_pretty_grid_name( self.g3d.name )
    ylim = [0.0, 0.0]
    # define axis labels and arrays
    if self.g3d.type == picdefs.hipace.h5.g3dtypes.density:
      ylim[0] = np.amin(self.line)
      ylim[1] = np.amax(self.line)
    elif self.g3d.type == picdefs.hipace.h5.g3dtypes.field:
      self.colormap = cm.coolwarm
      ylim[0] = -np.amax(abs(self.line))
      ylim[1] = np.amax(abs(self.line))
    elif self.g3d.type == picdefs.hipace.h5.g3dtypes.current:
      self.colormap = cm.coolwarm
      ylim[0] = -np.amax(abs(self.line))
      ylim[1] = np.amax(abs(self.line))
    self.ylim = ylim 

  def read( self ):
    # read line
    if self.args.lout_idx != None:
      lout_idx = list(self.args.lout_idx)

    if 'z' == self.args.loutax: 
      if self.args.lout_idx == None:
        # Default: central lineout
        idx1 = math.floor(self.g3d.nx[1]/2) - 1
        idx2 = math.floor(self.g3d.nx[2]/2) - 1
        self.line = self.g3d.read_line(i1=idx1, i2=idx2)
        if self.g3d.nx[1]%2 == 0 and self.g3d.nx[2]%2 == 0:
          line01 = self.g3d.read_line(i1=idx1, i2=idx2+1)
          line10 = self.g3d.read_line(i1=idx1+1, i2=idx2)
          line11 = self.g3d.read_line(i1=idx1+1, i2=idx2+1)
          self.line = ( self.line + line01 + line10 + line11 )/4
        elif self.g3d.nx[1]%2 == 1 and self.g3d.nx[2]%2 == 0:
          self.line = ( self.line + self.g3d.read_line(i1=idx1, i2=idx2+1) )/2
        elif self.g3d.nx[1]%2 == 0 and self.g3d.nx[2]%2 == 1:
          self.line = ( self.line + self.g3d.read_line(i1=idx1+1, i2=idx2) )/2
      else:
        self.line = self.g3d.read_slice(i1=lout_idx[0], i2=lout_idx[1])

    elif 'x' == self.args.loutax:
      if self.args.lout_idx == None:
        # Default: central lineout
        idx1 = math.floor(self.g3d.nx[0]/2) - 1
        idx2 = math.floor(self.g3d.nx[2]/2) - 1
        self.line = self.g3d.read_line(i0=idx1, i2=idx2)
        if self.g3d.nx[0]%2 == 0 and self.g3d.nx[2]%2 == 0:
          line01 = self.g3d.read_line(i0=idx1, i2=idx2+1)
          line10 = self.g3d.read_line(i0=idx1+1, i2=idx2)
          line11 = self.g3d.read_line(i0=idx1+1, i2=idx2+1)
          self.line = ( self.line + line01 + line10 + line11 )/4
        elif self.g3d.nx[0]%2 == 1 and self.g3d.nx[2]%2 == 0:
          self.line = ( self.line + self.g3d.read_line(i0=idx1, i2=idx2+1) )/2
        elif self.g3d.nx[0]%2 == 0 and self.g3d.nx[2]%2 == 1:
          self.line = ( self.line + self.g3d.read_line(i0=idx1+1, i2=idx2) )/2
      else:
        self.line = self.g3d.read_slice(i0=lout_idx[0], i2=lout_idx[1])

    elif 'y' == self.args.loutax:
      if self.args.lout_idx == None:
        # Default: central lineout
        idx1 = math.floor(self.g3d.nx[0]/2) - 1
        idx2 = math.floor(self.g3d.nx[1]/2) - 1
        self.line = self.g3d.read_line(i0=idx1, i1=idx2)
        if self.g3d.nx[0]%2 == 0 and self.g3d.nx[1]%2 == 0:
          line01 = self.g3d.read_line(i0=idx1, i1=idx2+1)
          line10 = self.g3d.read_line(i0=idx1+1, i1=idx2)
          line11 = self.g3d.read_line(i0=idx1+1, i1=idx2+1)
          self.line = ( self.line + line01 + line10 + line11 )/4
        elif self.g3d.nx[0]%2 == 1 and self.g3d.nx[1]%2 == 0:
          self.line = ( self.line + self.g3d.read_line(i0=idx1, i1=idx2+1) )/2
        elif self.g3d.nx[0]%2 == 0 and self.g3d.nx[1]%2 == 1:
          self.line = ( self.line + self.g3d.read_line(i0=idx1+1, i1=idx2) )/2
      else:
        self.line = self.g3d.read_slice(i0=lout_idx[0], i1=lout_idx[1])

  def plot( self, ifsave=True ):  
    if self.args.verbose: print('Generating line plot')    
    saveformat = self.args.file_format  
    filesuffix = '_%06.f' % (np.floor(self.g3d.time))
    
    if self.args.save_prefix != parsedefs.save_prefix.name:
      fileprefix = self.args.save_prefix
    else: 
      fileprefix = self.g3d.name

    savename = fileprefix + filesuffix + '.' + saveformat

    fig = plt.figure()
    cax = plt.plot( self.x_array, 
                    self.line)
    ax = plt.gca()
    ax.set_ylabel(self.ylabel, fontsize=14)
    ax.set_xlabel(self.xlabel, fontsize=14)

    self.mkdirs_if_nexist()

    if saveformat==parsedefs.file_format.png:
      fig.savefig(  self.args.savepath + '/' + savename, 
                format=saveformat,
                dpi=600)
    else:    
      fig.savefig(  self.args.savepath + '/' + savename, 
                    format=saveformat)
    if self.args.verbose: 
      print('Saved "' + savename + '" at: ' + self.args.savepath)    
    
    if self.args.ifshow: plt.show()
    plt.close(fig)



def plotfiles(args, ptype='none'):
  if not args.path:
    print('Error: No file provided!')
    sys.exit()

  for path in args.path:
    if os.path.isfile(path) :
      file = path
      if is_h5g3d_file(file): 
        if ptype == 'slice':
          g3d_p = G3d_plot_slice(file, args)
        elif ptype == 'line':
          g3d_p = G3d_plot_line(file, args)
        g3d_p.plot()
      else:
        print('Skipping: ' + file)  
    elif os.path.isdir(path):
      print('"' + path + '"' + ' is a directory.')
      if args.process_all == True:
        print('Processing all g3d files in the provided directory.')
        for root, dirs, files in os.walk(path):  
          for filename in files:
            file = root + '/' + filename
            if is_h5g3d_file(file): 
              if ptype == 'slice':
                g3d_p = G3d_plot_slice(file, args)
              elif ptype == 'line':
                g3d_p = G3d_plot_line(file, args)
              g3d_p.plot()
            else:
              print('Skipping: ' + file)        
        sys.exit()
      else:
        print('Error: Use the flag "-a" to process all files in the provided directory!')
        sys.exit() 
    elif not os.path.exists(path):
      print('Error: Provided path does not exist!')
      sys.exit()    
    else:
      print('Error: Provided path is neither a file nor a directory!')
      sys.exit()
