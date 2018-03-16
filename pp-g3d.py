#!/usr/bin/env python3

import sys
import os
import math
import argparse
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import MaxNLocator
from matplotlib import cm
import pp_defs
from pp_h5dat import Grid3d
from pp_h5dat import H5FList
from pp_h5dat import mkdirs_if_nexist

# Parse defaults/definitions
class parsedefs:
    class file_format:
        png = 'png'
        eps = 'eps'
        pdf = 'pdf'
    class zax:
        zeta = 'zeta'
        z = 'z'
        xi = 'xi'
    class save:
        prefix = 'g3d_name'
        path = './plots'

def g3d_parser():

    desc = """This is the picpy postprocessing tool."""

    parser = argparse.ArgumentParser( add_help=False,
                                      description=desc)
    parser.add_argument(  'path',
                          metavar = 'PATH',
                          nargs = '*',
                          help = 'Path to grid file.')
    parser.add_argument(  "-v", "--verbose",
                          dest = "verbose",
                          action="store_true",
                          default=True,
                          help = "Print info (default: %(default)s).")
    parser.add_argument(  "-q", "--quiet",
                          dest = "verbose",
                          action = "store_false",
                          help = "Don't print info.")
    parser.add_argument(  "--show",
                          dest = "ifshow",
                          action = "store_true",
                          default = False,
                          help = "Show figure (default: %(default)s).")
    parser.add_argument(  "--name-prefix",
                          action = "store",
                          dest = "save_prefix",
                          metavar = "NAME",
                          default = parsedefs.save.prefix,
                          help = """Define customized prefix of output filename (default: %(default)s).""")
    parser.add_argument(  "-c", "--code",
                          action="store",
                          dest="piccode",
                          choices = [pp_defs.code.hipace, pp_defs.code.osiris,],
                          default = pp_defs.code.hipace,
                          help="PIC code (default: %(default)s).")
    parser.add_argument(  "-z", "--z-axis",
                          action='store',
                          dest="zax",
                          choices=[ parsedefs.zax.zeta,
                                  parsedefs.zax.z,
                                  parsedefs.zax.xi,],
                          default=parsedefs.zax.zeta,
                          help= "z-axis type (default: %(default)s).")
    return parser

def g3d_slice_subparser(subparsers, parent_parser):
    parser = subparsers.add_parser( "slice", parents=[parent_parser],
                                    help="Grid 3D slice plotting")    
    # Slice plot specific arguments
    parser.add_argument(  "-p", "--plane",
                          action='store',
                          dest="plane",
                          choices=[ 'xy', 'yz', 'xz',
                                    'yx', 'zy', 'zx'],
                          default='zx',
                          help= """Plane to be plotted (default: %(default)s).""")
    parser.add_argument(  "--pidx",
                          action='store',
                          dest="plane_index",
                          metavar="PLANE-IDX",
                          default=None,
                          type=int,
                          help= """Index of plane.""")
    parser.add_argument(  "--ppos",
                          action='store',
                          dest="plane_pos",
                          metavar="PLANE-POS",
                          default=None,
                          type=float,
                          help= """Position of plane.""")    
    parser.add_argument(  "--cscale",
                          default="lin",
                          dest="cscale",
                          choices=[ "lin", "log",],
                          help= "z-axis type (default: %(default)s).")
    parser.add_argument(  '--cblim',
                          help='Colorbar axis limits',
                          action='store',
                          dest="cblim",
                          metavar=('CBMIN', 'CBMAX'),
                          nargs=2,
                          type=float,
                          default=None)
    parser.add_argument(  "--save-path",
                          action="store",
                          dest="savepath",
                          metavar="PATH",
                          default=parsedefs.save.path + '/g3d-slice',
                          help = """Path to which generated files will be saved.
                              (default: %(default)s)""")
    parser.add_argument(  "-f", "--format",
                          action='store',
                          dest="file_format",
                          choices=[ 'png',
                                    'pdf',
                                    'eps',],
                          default='png',
                          help= """Format of output file (default: %(default)s).""")    
    return parser


def g3d_line_subparser(subparsers, parent_parser):
    parser = subparsers.add_parser( "line", parents=[parent_parser],
                                    help="Grid 3D line plotting")    
    # Line plot specific arguments
    parser.add_argument(  "-a", "--axis",
                          action='store',
                          dest="lineax",
                          choices=[ 'x', 'y', 'z'],
                          default='z',
                          help= """Axis along which lineout is generated (default: %(default)s).""")
    parser.add_argument(  '-i', '--idx',
                          help='Indices for which lineout is taken.',
                          action='store',
                          dest="lout_idx",
                          metavar=('idx0', 'idx1'),
                          nargs=2,
                          type=int,
                          default=None)
    parser.add_argument(  "--zeta-pos",
                          action='store',
                          dest="lout_zeta_pos",
                          metavar="LOUT-ZETA-POS",
                          type=float,    
                          default=None,
                          help= """Zeta-position at which lineout is generated (default: %(default)s).""")    
    parser.add_argument(  "-s", "--save-path",
                          action="store",
                          dest="savepath",
                          metavar="PATH",
                          default=parsedefs.save.path + '/g3d-line',
                          help = """Path to which generated files will be saved.
                              (default: %(default)s)""")
    parser.add_argument(  "-f", "--format",
                          action='store',
                          dest="file_format",
                          choices=[ 'png',
                                    'pdf',
                                    'eps',],
                          default='eps',
                          help= """Format of output file (default: %(default)s).""")
    parser.add_argument(  '-r','--range',
                          help='Range of lineout.',
                          action='store',
                          dest="range",
                          metavar=('XMIN', 'XMAX'),
                          nargs=2,
                          type=float,
                          default=None)
    parser.add_argument(  "--lino",
                          dest = "if_lineout",
                          action="store_true",
                          default=True,
                          help = "Generate lineout (default: %(default)s).")
    parser.add_argument(  "--int",
                          dest = "if_lineout",
                          action = "store_false",
                          help = "Generate integral (default: False).")                         
    return parser

# Converting HDF strings of grid quantity namnes
# to LaTeX strings
def gen_pretty_grid_name( gname ):
    if gname == 'ExmBy':
        return r'$(E_x-B_y)/E_0$'
    elif gname == 'EypBx':
        return r'$(E_y+B_x)/E_0$'
    elif gname == 'Ez':
        return r'$E_z/E_0$'
    elif gname == 'Bx':
        return r'$B_x/E_0$'
    elif gname == 'By':
        return r'$B_y/E_0$'
    elif gname == 'Bz':
        return r'$B_z/E_0$'
    elif gname == 'Jx':
        return r'$J_x/cn_0$'
    elif gname == 'Jy':
        return r'$J_y/cn_0$'
    elif gname == 'Jz':
        return r'$J_z/cn_0$'
    elif gname == 'Jz':
        return r'$J_z$'
    elif gname == 'mpsi':
        return r'$-\psi/m_e c^2$'
    elif 'plasma' in gname:
        if 'charge' in gname:
            return r'$\rho_p/e n_0$'
        else:
            return r'$n_p/n_0$'
    elif 'beam' in gname:
        if 'charge' in gname:
            return r'$\rho_b/e n_0$'
        else:
            return r'$n_b/n_0$'                   
    else:
        return gname


# Returning boolean: if file extension is hdf5 extension
def is_h5_file(fext):
    return any(fext == h5ext for h5ext in pp_defs.fexts.hdf5)

# Returning boolean: if file name contains name of grid quantity
# ['density', 'field', 'current']
def is_g3d_file(fname):
    return any((mq in fname) for mq in pp_defs.hipace.h5.g3dtypes.list)

# Returning boolean:  if file extension is hdf5 extension and
#                     if file name contains name of grid quantity
def is_h5g3d_file(file):
    fname, fext = os.path.splitext(file)
    return is_h5_file(fext) and is_g3d_file(fname)

class G3d_hslab:
    def __init__(self, file, i0=None, i1=None, i2=None):
        self.file = file
        # Do everything here that does not require any args!!!
        # - readin of hyperslab
        # - axes arrays
        # etc.       

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

    def is_number_density( self ):
        name = self.g3d.name
        if 'plasma' in name:
            if 'charge' in name:
                return False
            else:
                return True
        elif 'beam' in name:
            if 'charge' in name:
                return False
            else:
                return True                  
        else:
            return False

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
                if (self.args.plane_index == None) and (self.args.plane_pos == None):
                    index = math.floor(self.g3d.nx[2]/2) - 1
                    self.slice = self.g3d.read(i2=index)
                    if self.g3d.nx[2]%2 == 0:
                        self.slice = (self.slice + self.g3d.read(i2=index+1))/2
                elif (self.args.plane_index != None) and (self.args.plane_pos == None):
                    self.slice = self.g3d.read(i2=self.args.plane_index)
                elif (self.args.plane_index == None) and (self.args.plane_pos != None): 
                    self.slice = self.g3d.read(x2=self.args.plane_pos)
                else:
                    print('ERROR: plane-index can''t be used in conjunction with plane-position!')
                    sys.exit(1) 

            elif 'y' in self.args.plane:
                if (self.args.plane_index == None) and (self.args.plane_pos == None):
                    index = math.floor(self.g3d.nx[1]/2) - 1
                    self.slice = self.g3d.read(i1=index)
                    if self.g3d.nx[1]%2 == 0:
                        self.slice = ( self.slice + self.g3d.read(i1=index+1) )/2
                elif (self.args.plane_index != None) and (self.args.plane_pos == None): 
                    self.slice = self.g3d.read(i1=self.args.plane_index)
                elif (self.args.plane_index == None) and (self.args.plane_pos != None): 
                    self.slice = self.g3d.read(x1=self.args.plane_pos)                    
                else:
                    print('ERROR: plane-index can''t be used in conjunction with plane-position!')
                    sys.exit(1) 
        elif ('x' in self.args.plane) and ('y' in self.args.plane):
            if (self.args.plane_index == None) and (self.args.plane_pos == None):
                index = math.floor(self.g3d.nx[0]/2) - 1
                self.slice = self.g3d.read(i0=index)
                if self.g3d.nx[0]%2 == 0:
                    self.slice = ( self.slice + self.g3d.read(i0=index+1) )/2
            elif (self.args.plane_index != None) and (self.args.plane_pos == None):
                self.slice = self.g3d.read(i0=self.args.plane_index)
            elif (self.args.plane_index == None) and (self.args.plane_pos != None):    
                self.slice = self.g3d.read(x0=self.args.plane_pos) 
            else:
                print('ERROR: plane-index can''t be used in conjunction with plane-position!')
                sys.exit(1)                 

        if self.args.plane in ['xy','zx','zy']:
            self.slice = np.transpose( self.slice )

        if self.is_number_density():
            self.slice = np.abs(self.slice)

    def set_cmap( self ):
        if self.args.cscale == "log":
            self.slice = np.log(abs(self.slice))

        cblim = [0.0, 0.0]

        if self.g3d.type == pp_defs.hipace.h5.g3dtypes.density:
            if self.is_number_density():
                self.colormap = cm.PuBu;
                cblim[0] = np.amin(self.slice)
                cblim[1] = np.amax(self.slice)
                if (cblim[0] == 0.0) and (cblim[1] == 0.0):
                    cblim[0] = 0.0
                    cblim[1] = 1.0
            else:
                self.colormap = cm.RdGy;
                cblim[0] = -np.amax(abs(self.slice))
                cblim[1] = np.amax(abs(self.slice))              
                if (cblim[0] == 0.0) and (cblim[1] == 0.0):
                    cblim[0] = -1.0
                    cblim[1] = 1.0

        elif self.g3d.type == pp_defs.hipace.h5.g3dtypes.field:
            self.colormap = cm.RdBu
            cblim[0] = -np.amax(abs(self.slice))
            cblim[1] = np.amax(abs(self.slice))
            if (cblim[0] == 0.0) and (cblim[1] == 0.0):
                cblim[0] = -1.0
                cblim[1] = 1.0            
        elif self.g3d.type == pp_defs.hipace.h5.g3dtypes.current:
            self.colormap = cm.PuOr
            cblim[0] = -np.amax(abs(self.slice))
            cblim[1] = np.amax(abs(self.slice))
            if (cblim[0] == 0.0) and (cblim[1] == 0.0):
                cblim[0] = -1.0
                cblim[1] = 1.0 

        if self.args.cblim != None:
            cblim = list(self.args.cblim)

        self.cblim = cblim

    def plot( self, ifsave=True ):
        if self.args.verbose: print('Generating slice plot')
        saveformat = self.args.file_format
        filesuffix = '_%06.f' % (np.floor(self.g3d.time))

        if self.args.save_prefix != parsedefs.save.prefix:
            fileprefix = self.args.save_prefix
        else:
            fileprefix = self.g3d.name

        if self.g3d.is_subgrid():
            sg_str = '_subgrid'
        else:
            sg_str = ''   

        savename = "%s%s_%s%s.%s" % (fileprefix, sg_str, self.args.plane, filesuffix, saveformat)

        levels = MaxNLocator(nbins=256).tick_values(self.cblim[0], self.cblim[1])

        fig = plt.figure()
        cax = plt.contourf( self.x_array,
                            self.y_array,
                            self.slice,
                            levels=levels,
                            vmin=self.cblim[0], vmax=self.cblim[1],
                            cmap=self.colormap)
        ax = plt.gca()
        ax.set_ylabel(self.ylabel, fontsize=14)
        ax.set_xlabel(self.xlabel, fontsize=14)
        cbar = fig.colorbar(cax)
        cbar.ax.set_ylabel( gen_pretty_grid_name( self.g3d.name ), fontsize=14 )

        axpos1 = ax.get_position() # get the original position 
        cbaxpos1 = cbar.ax.get_position() # get the original position 
        xshift = 0.05;
        wcompress = 1.05
        axpos2 = [axpos1.x0 + xshift, axpos1.y0,  (axpos1.width - xshift) / wcompress, axpos1.height ]
        ax.set_position(axpos2) # set a new position 
        axpos2 = ax.get_position()
        cbaxpos2 = [cbaxpos1.x0 + (axpos2.x0 - axpos1.x0) + (axpos2.width - axpos1.width), \
                    cbaxpos1.y0, \
                    cbaxpos1.width, \
                    cbaxpos1.height ] 
  
        cbar.ax.set_position(cbaxpos2)

        mkdirs_if_nexist(self.args.savepath)

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

        self.set_xaxis( self.args.lineax )
        if self.args.verbose: print('Reading data...')
        self.read()
        if self.args.verbose: print('Read-in completed.')
        self.set_yaxis()

    def set_yaxis( self ):
        if self.args.if_lineout:
            self.ylabel = gen_pretty_grid_name( self.g3d.name )
        else:
            self.ylabel = r'$k_p^2 \int \int$' + gen_pretty_grid_name( self.g3d.name )
            if self.args.lineax == 'z':
                self.ylabel = self.ylabel + r'$\,dx dy$'
            if self.args.lineax == 'x':
                self.ylabel = self.ylabel + r'$\,dy dz$'
            if self.args.lineax == 'y':
                self.ylabel = self.ylabel + r'$\,dx dz$'
        ylim = [0.0, 0.0]
        # define axis labels and arrays
        if self.g3d.type == pp_defs.hipace.h5.g3dtypes.density:
            ylim[0] = np.amin(self.line)
            ylim[1] = np.amax(self.line)
        elif self.g3d.type == pp_defs.hipace.h5.g3dtypes.field:
            self.colormap = cm.coolwarm
            ylim[0] = -np.amax(abs(self.line))
            ylim[1] = np.amax(abs(self.line))
        elif self.g3d.type == pp_defs.hipace.h5.g3dtypes.current:
            self.colormap = cm.coolwarm
            ylim[0] = -np.amax(abs(self.line))
            ylim[1] = np.amax(abs(self.line))
        self.ylim = ylim

    def read( self ):
        # read line
        if self.args.lout_idx != None:
            lout_idx = list(self.args.lout_idx)

        if 'z' == self.args.lineax:
            if self.args.if_lineout:
                if self.args.lout_idx == None:
                    # Default: central lineout
                    idx1 = math.floor(self.g3d.nx[1]/2) - 1
                    idx2 = math.floor(self.g3d.nx[2]/2) - 1
                    self.line = self.g3d.read(i1=idx1, i2=idx2)
                    if self.g3d.nx[1]%2 == 0 and self.g3d.nx[2]%2 == 0:
                        line01 = self.g3d.read(i1=idx1, i2=idx2+1)
                        line10 = self.g3d.read(i1=idx1+1, i2=idx2)
                        line11 = self.g3d.read(i1=idx1+1, i2=idx2+1)
                        self.line = ( self.line + line01 + line10 + line11 )/4
                    elif self.g3d.nx[1]%2 == 1 and self.g3d.nx[2]%2 == 0:
                        self.line = ( self.line + self.g3d.read(i1=idx1, i2=idx2+1) )/2
                    elif self.g3d.nx[1]%2 == 0 and self.g3d.nx[2]%2 == 1:
                        self.line = ( self.line + self.g3d.read(i1=idx1+1, i2=idx2) )/2
                else:
                    self.line = self.g3d.read(i1=lout_idx[0], i2=lout_idx[1])
            else:
                self.line = self.g3d.read_integrate(ax1=True,ax2=True)

        elif 'x' == self.args.lineax:
            if self.args.if_lineout:
                if (self.args.lout_idx == None) and (self.args.lout_zeta_pos == None):
                    # Default: central lineout
                    idx1 = math.floor(self.g3d.nx[0]/2) - 1
                    idx2 = math.floor(self.g3d.nx[2]/2) - 1
                    self.line = self.g3d.read(i0=idx1, i2=idx2)
                    if self.g3d.nx[0]%2 == 0 and self.g3d.nx[2]%2 == 0:
                        line01 = self.g3d.read(i0=idx1, i2=idx2+1)
                        line10 = self.g3d.read(i0=idx1+1, i2=idx2)
                        line11 = self.g3d.read(i0=idx1+1, i2=idx2+1)
                        self.line = ( self.line + line01 + line10 + line11 )/4
                    elif self.g3d.nx[0]%2 == 1 and self.g3d.nx[2]%2 == 0:
                        self.line = ( self.line + self.g3d.read(i0=idx1, i2=idx2+1) )/2
                    elif self.g3d.nx[0]%2 == 0 and self.g3d.nx[2]%2 == 1:
                        self.line = ( self.line + self.g3d.read(i0=idx1+1, i2=idx2) )/2
                elif (self.args.lout_idx != None) and (self.args.lout_zeta_pos == None):
                    self.line = self.g3d.read(i0=lout_idx[0], i2=lout_idx[1])
                elif (self.args.lout_idx == None) and (self.args.lout_zeta_pos != None):
                    self.line = self.g3d.read(x0=self.args.lout_zeta_pos, x2=0.0)
                else:
                    print('ERROR: lineout-index can''t be used in conjunction with lineout-zeta-position!')
                    sys.exit(1)
            else:
                self.line = self.g3d.read_integrate(ax0=True,ax2=True)


        elif 'y' == self.args.lineax:
            if self.args.if_lineout:
                if (self.args.lout_idx == None) and (self.args.lout_zeta_pos == None):
                    # Default: central lineout
                    idx1 = math.floor(self.g3d.nx[0]/2) - 1
                    idx2 = math.floor(self.g3d.nx[1]/2) - 1
                    self.line = self.g3d.read(i0=idx1, i1=idx2)
                    if self.g3d.nx[0]%2 == 0 and self.g3d.nx[1]%2 == 0:
                        line01 = self.g3d.read(i0=idx1, i1=idx2+1)
                        line10 = self.g3d.read(i0=idx1+1, i1=idx2)
                        line11 = self.g3d.read(i0=idx1+1, i1=idx2+1)
                        self.line = ( self.line + line01 + line10 + line11 )/4
                    elif self.g3d.nx[0]%2 == 1 and self.g3d.nx[1]%2 == 0:
                        self.line = ( self.line + self.g3d.read(i0=idx1, i1=idx2+1) )/2
                    elif self.g3d.nx[0]%2 == 0 and self.g3d.nx[1]%2 == 1:
                        self.line = ( self.line + self.g3d.read(i0=idx1+1, i1=idx2) )/2
                elif (self.args.lout_idx != None) and (self.args.lout_zeta_pos == None):
                    self.line = self.g3d.read(i0=lout_idx[0], i1=lout_idx[1])
                elif (self.args.lout_idx == None) and (self.args.lout_zeta_pos != None):
                    self.line = self.g3d.read(x0=self.args.lout_zeta_pos, x1=0.0)
                else:
                    print('ERROR: lineout-index can''t be used in conjunction with lineout-zeta-position!')
                    sys.exit(1)
            else:
                self.line = self.g3d.read_integrate(ax0=True,ax1=True)                                

        if self.is_number_density():
            self.line = np.abs(self.line)

    def plot( self, ifsave=True ):
        if self.args.verbose: print('Generating line plot')
        saveformat = self.args.file_format
        filesuffix = '_%06.f' % (np.floor(self.g3d.time))

        if self.args.save_prefix != parsedefs.save.prefix:
            fileprefix = self.args.save_prefix
        else:
            fileprefix = self.g3d.name

        if self.g3d.is_subgrid():
            sg_str = '_subgrid'
        else:
            sg_str = ''   

        if self.args.if_lineout:
            integral_str = ''
        else:
            integral_str = '_int'  

        savename = "%s%s_%s%s%s.%s" % (fileprefix, \
                                       sg_str, \
                                       self.args.lineax, \
                                       integral_str, \
                                       filesuffix, \
                                       saveformat)

        fig = plt.figure()
        cax = plt.plot( self.x_array,
                        self.line)
        ax = plt.gca()
        ax.set_ylabel(self.ylabel, fontsize=14)
        ax.set_xlabel(self.xlabel, fontsize=14)

        max_abs_val = np.max(np.abs(self.line))

        if max_abs_val > 0.0:
            if not (-3.0 < math.log(max_abs_val,10) < 3.0):
                ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
                plt.gcf().subplots_adjust(left=0.18)
            else:
                plt.gcf().subplots_adjust(left=0.15)      

        if self.args.range != None:
            ax.set_xlim(self.args.range)
            lims = ax.get_xlim()
            idx = np.where( (self.x_array >= lims[0]) & (self.x_array <= lims[1]) )[0]
            ax.set_ylim( self.line[idx].min(), self.line[idx].max() )

        mkdirs_if_nexist(self.args.savepath)

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


def slice(args):

    h5flist = H5FList(args.path, h5ftype='g3d')
    flist = h5flist.get()

    for file in flist:
        g3d_p = G3d_plot_slice(file, args)
        g3d_p.plot()        


def line(args):

    h5flist = H5FList(args.path, h5ftype='g3d')
    flist = h5flist.get()

    for file in flist:
        g3d_p = G3d_plot_line(file, args)
        g3d_p.plot()  



def main():
    parser = argparse.ArgumentParser()
    g3d_subparsers = parser.add_subparsers(title="plot-type")

    g3d_parent_parser = g3d_parser()

    g3d_ssp = g3d_slice_subparser(  subparsers=g3d_subparsers,
                                    parent_parser=g3d_parent_parser)
    g3d_ssp.set_defaults(func=slice)

    g3d_lsp = g3d_line_subparser( subparsers=g3d_subparsers,
                                  parent_parser=g3d_parent_parser)
    g3d_lsp.set_defaults(func=line)

    args = parser.parse_args()
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    args.func(args)


if __name__ == "__main__":
    main()         