#!/usr/bin/env python3

import sys
import os
import math
import argparse
import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import MaxNLocator
from matplotlib import cm

mypath = os.path.dirname(os.path.abspath( __file__ ))
incpath = os.path.split(mypath)[0] + '/inc'
sys.path.append(incpath)
import pp_defs
from pp_h5dat import Grid3d
from pp_h5dat import Grid2d
from pp_h5dat import H5FList
from pp_h5dat import H5Plot
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
                          choices = [pp_defs.code.hipace, pp_defs.code.osiris, pp_defs.code.inferno],
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
    parser.add_argument(  "--int",
                          dest = "if_integrate",
                          action = "store_true",
                          default=False,                          
                          help = "Plot integrated quantity (default: %(default)s).")
    parser.add_argument(  "--smooth",
                          dest = "smoothing",
                          action = "store",
                          default = False,
                          type=float,
                          metavar="SIGMA",
                          help = "Smoothing of the line plot. Sigma is the input "
                          " into the 1D gaussian filter (default: %(default)s).")
    parser.add_argument(  "--avgax",
                          action='store',
                          dest="avgax",
                          choices=[ 'x', 'y', 'z'],
                          default=None,
                          help= """Axis-position averaged for (default: %(default)s).""")    
    parser.add_argument(  "--xlim",
                          help='Customize x-axis limits',
                          action='store',
                          dest="xlim",
                          metavar=('xmin', 'xmax'),
                          type=float,
                          nargs=2,
                          default=None)  
    parser.add_argument(  "--ylim",
                          help='Customize y-axis limits',
                          action='store',
                          dest="ylim",
                          metavar=('ymin', 'ymax'),
                          type=float,
                          nargs=2,
                          default=None)
    parser.add_argument(  "--ylabelpad",
                        help='Customize the label pad for the second y-axis',
                        action='store',
                        dest="ylabelpad",
                        type=int, 
                        default=None)   
    parser.add_argument(  "--h5off",
                          dest = "h5plot_off",
                          action="store_true",
                          default=False,
                          help = "Save plot as hdf5 file (Default: %(default)s).")
    parser.add_argument(  "-s", "--save-path",
                          action="store",
                          dest="savepath",
                          metavar="PATH",
                          default=None,
                          help = """Path to which generated files will be saved.
                              (default: %(default)s)""")
    parser.add_argument(  "--latexfont",
                          dest = "latexfont",
                          action="store_true",
                          default=False,
                          help = "Use LaTeX font (Default: %(default)s).")
    parser.add_argument(  "--fontsize",
                          dest = "fontsize",
                          action="store",
                          default=14,
                          type=int,
                          help = "Set fontsize of axis labels (Default: %(default)s).")
    parser.add_argument(  "--cbarpad",
                          dest = "cbarpad",
                          action="store",
                          default=0,
                          type=int,
                          help = "Set distance of title from cbar (Default: %(default)s).")
    parser.add_argument(  "--axticklabelsize",
                        dest = "axticklabelsize",
                        action="store",
                        default=11,
                        type=int,
                        help = "Set fontsize of axis tick labels (Default: %(default)s).")
    parser.add_argument(  "--dpi",
                          action='store',
                          dest="dpi",
                          default=300,
                          type=int,
                          help= """Dots per inch for png output (default: %(default)s).""")  
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
    parser.add_argument(  "--clog",
                          dest = "clog",
                          action = "store_true",
                          default = False,
                          help = "Log color scale (default: %(default)s).")
    parser.add_argument(  '--clim',
                          help='Colorbar axis limits',
                          action='store',
                          dest="clim",
                          metavar=('CBMIN', 'CBMAX'),
                          nargs=2,
                          type=float,
                          default=None)  
    parser.add_argument(  "-f", "--format",
                          action='store',
                          dest="file_format",
                          choices=[ 'png',
                                    'pdf',
                                    'eps',],
                          default='png',
                          help= """Format of output file (default: %(default)s).""")
    parser.add_argument(  "--manual",
                          dest = "manual",
                          action="store_true",
                          default=False,
                          help = "Enables manual multiplots (timestamps don't have to agree) (Default: %(default)s).")
    parser.add_argument(  "--data2",
                        dest = "data2",
                        nargs = '*',
                        default = False,
                        help = "Path to the second data set which shall be overlaid over the first one. \n")
    parser.add_argument(  "--clog2",
                          dest = "clog2",
                          action = "store_true",
                          default = False,
                          help = "Log color scale for the second data set (default: %(default)s).")
    parser.add_argument(  "--clim2",
                        help='Colorbar axis limits for the second data set',
                        action='store',
                        dest="clim2",
                        metavar=('CBMIN', 'CBMAX'),
                        nargs=2,
                        type=float,
                        default=None)
    parser.add_argument(  "--cm2",
                          dest = "cm2",
                          action = "store",
                          default = False,
                          help = "Choose different than standard matplotlib colormap for the second data set. Example: YlGn (default: %(default)s).")
    parser.add_argument(  "--data3",
                        dest = "data3",
                        nargs = '*',
                        default = False,
                        help = "Path to the third data set which shall be overlaid over the first one. \n")
    parser.add_argument(  "--clog3",
                          dest = "clog3",
                          action = "store_true",
                          default = False,
                          help = "Log color scale for the third data set (default: %(default)s).")
    parser.add_argument(  "--diff",
                          dest = "diff",
                          action = "store_true",
                          default = False,
                          help = "Plots the difference between the first and the second data set (default: %(default)s).")
    parser.add_argument(  "--clim3",
                        help='Colorbar axis limits for the third data set',
                        action='store',
                        dest="clim3",
                        metavar=('CBMIN', 'CBMAX'),
                        nargs=2,
                        type=float,
                        default=None)
    parser.add_argument(  "--cm3",
                          dest = "cm3",
                          action = "store",
                          default = False,
                          help = "Choose different than standard matplotlib colormap for the third data set. Example: YlGn (default: %(default)s).")
    parser.add_argument(  "--ptype",
                          default="pcolormesh",
                          dest="ptype",
                          choices=[ "pcolor", "pcolormesh", "imshow", "pcolorfast", "contourf"],
                          help= "Plot color type (default: %(default)s).")
    parser.add_argument(  "--line_data",
                        dest = "line_data",
                        nargs = '*',
                        default = False,
                        help = "Path to the line data set, which shall be overlaid over the 2D plot. \n")
    parser.add_argument(  "-a", "--axis",
                          action='store',
                          dest="lineax",
                          choices=[ 'x', 'y', 'z'],
                          default='z',
                          help= """Axis along which lineout is generated for the line data overlay plot (default: %(default)s).""")
    parser.add_argument(  '-i', '--idx',
                          help='Indices for which lineout is taken in the line plot.',
                          action='store',
                          dest="lout_idx",
                          metavar=('IDX0', 'IDX1'),
                          nargs=2,
                          type=int,
                          default=None)
    parser.add_argument(  "--ylog",
                          dest = "absylog",
                          action="store_true",
                          default=False,
                          help = "Plot abs log of y-data in line out overlay (default: %(default)s).")
    parser.add_argument(  '-r','--range',
                          help='Range of lineout.',
                          action='store',
                          dest="range",
                          metavar=('XMIN', 'XMAX'),
                          nargs=2,
                          type=float,
                          default=None)
    parser.add_argument(  "--ylim2",
                        help='Customize y-axis limits for the line plot',
                        action='store',
                        dest="ylim2",
                        metavar=('ymin', 'ymax'),
                        type=float,
                        nargs=2,
                        default=False)
    parser.add_argument(  "--y2labelpad",
                        help='Customize the label pad for the second y-axis',
                        action='store',
                        dest="y2labelpad",
                        type=int, 
                        default=None)       
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
                          metavar=('IDX0', 'IDX1'),
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
    parser.add_argument(  '-r','--range',
                          help='Range of lineout.',
                          action='store',
                          dest="range",
                          metavar=('XMIN', 'XMAX'),
                          nargs=2,
                          type=float,
                          default=None) 
    parser.add_argument(  "--ylog",
                          dest = "absylog",
                          action="store_true",
                          default=False,
                          help = "Plot abs log of y-data (default: %(default)s).")                           
    parser.add_argument(  "-f", "--format",
                          action='store',
                          dest="file_format",
                          choices=[ 'png',
                                    'pdf',
                                    'eps',],
                          default='pdf',
                          help= """Format of output file (default: %(default)s).""")
    parser.add_argument(  "--ylim2",
                        help='Customize y-axis limits for the line plot'
                                'Use only for overlay plots!',
                        action='store',
                        dest="ylim2",
                        metavar=('ymin', 'ymax'),
                        type=float,
                        nargs=2,
                        default=False)   
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

def round_figures(x, n): 
    """Returns x rounded to n significant figures."""
    return round(x, int(n - math.ceil(math.log10(abs(x)))))

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
        self.app_str = ''
        self.root_savepath = './plots'

        # Reading hdf5 attributes:
        if self.args.verbose:  
            sys.stdout.write('Getting attributes of %s\n' % file)
            sys.stdout.flush()

        if self.args.piccode == pp_defs.code.inferno:
            self.g3d = Grid2d(file)
        else:
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
                sys.stdout.write('Error: No/wrong z-axis option selected!\n')
                sys.stdout.flush()
                sys.exit()
        elif xax_str == 'x':
            self.x_array = self.g3d.get_x_arr(1)
            self.xlabel = r'$k_p x$'
        elif xax_str == 'y':
            self.x_array = self.g3d.get_x_arr(2)
            self.xlabel = r'$k_p y$'
        else:
            sys.stdout.write('Error: Wrong x-axis string!\n')
            sys.stdout.flush()
            sys.exit()

    def is_number_density( self ):
        name = self.g3d.get_name()
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

        self.__relsavepath = '/g3d-slice'

        self.set_xaxis( self.args.plane[0] )
        self.set_yaxis( self.args.plane[1] )
        if self.args.verbose: 
            sys.stdout.write('Reading data...\n')
            sys.stdout.flush()
        self.read()
        if self.args.verbose: 
            sys.stdout.write('Read-in completed.\n')
            sys.stdout.flush()
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
                sys.stdout.write('Error: No/wrong z-axis option selected!\n')
                sys.stdout.flush()
                sys.exit()
        elif yax_str == 'x':
            self.y_array = self.g3d.get_x_arr(1)
            self.ylabel = r'$k_p x$'
        elif yax_str == 'y':
            self.y_array = self.g3d.get_x_arr(2)
            self.ylabel = r'$k_p y$'
        else:
            sys.stdout.write('Error: Wrong y-axis string!\n')
            sys.stdout.flush()
            sys.exit()

    def read( self ):
        # read slice
        if 'z' in self.args.plane:
            if 'x' in self.args.plane:
                if self.args.if_integrate:
                    self.slice = self.g3d.read_integrate(ax2=True)
                else:
                    if self.args.piccode == pp_defs.code.inferno:
                        self.slice = self.g3d.read()
                    else:
                        if (self.args.plane_index == None) and (self.args.plane_pos == None):
                            index = math.floor(self.g3d.get_nx(2)/2) - 1
                            self.slice = self.g3d.read(i2=index)
                            if self.g3d.get_nx(2)%2 == 0:
                                self.slice = (self.slice + self.g3d.read(i2=index+1))/2
                        elif (self.args.plane_index != None) and (self.args.plane_pos == None):
                            self.slice = self.g3d.read(i2=self.args.plane_index)
                        elif (self.args.plane_index == None) and (self.args.plane_pos != None): 
                            self.slice = self.g3d.read(x2=self.args.plane_pos)
                        else:
                            sys.stdout.write('ERROR: plane-index can''t be used in conjunction with plane-position!\n')
                            sys.stdout.flush()
                            sys.exit(1)                    

            elif 'y' in self.args.plane:
                if self.args.if_integrate:
                    self.slice = self.g3d.read_integrate(ax1=True)
                else:    
                    if (self.args.plane_index == None) and (self.args.plane_pos == None):
                        index = math.floor(self.g3d.get_nx(1)/2) - 1
                        self.slice = self.g3d.read(i1=index)
                        if self.g3d.get_nx(1)%2 == 0:
                            self.slice = ( self.slice + self.g3d.read(i1=index+1) )/2
                    elif (self.args.plane_index != None) and (self.args.plane_pos == None): 
                        self.slice = self.g3d.read(i1=self.args.plane_index)
                    elif (self.args.plane_index == None) and (self.args.plane_pos != None): 
                        self.slice = self.g3d.read(x1=self.args.plane_pos)                    
                    else:
                        sys.stdout.write('ERROR: plane-index can''t be used in conjunction with plane-position!\n')
                        sys.stdout.flush()
                        sys.exit(1)                    

        elif ('x' in self.args.plane) and ('y' in self.args.plane):
            if self.args.if_integrate:
                self.slice = self.g3d.read_integrate(ax0=True)
            else:
                if (self.args.plane_index == None) and (self.args.plane_pos == None):
                    index = math.floor(self.g3d.get_nx(0)/2) - 1
                    self.slice = self.g3d.read(i0=index)
                    if self.g3d.get_nx(0)%2 == 0:
                        self.slice = ( self.slice + self.g3d.read(i0=index+1) )/2
                elif (self.args.plane_index != None) and (self.args.plane_pos == None):
                    self.slice = self.g3d.read(i0=self.args.plane_index)
                elif (self.args.plane_index == None) and (self.args.plane_pos != None):    
                    self.slice = self.g3d.read(x0=self.args.plane_pos) 
                else:
                    sys.stdout.write('ERROR: plane-index can''t be used in conjunction with plane-position!\n')
                    sys.stdout.flush()
                    sys.exit(1)                               

        if self.args.plane in ['xy','zx','zy']:
            self.slice = np.transpose( self.slice )

        if self.is_number_density():
            self.slice = np.abs(self.slice)

    def set_cmap( self, update=False, clim_input=None, cmap_input=False ):

        clim = [0.0, 0.0]
        if self.g3d.get_type() == pp_defs.hipace.h5.g3dtypes.density:
            if self.is_number_density():
                if 'beam' in self.g3d.get_name():
                    if not update:
                        self.colormap = cm.Reds;
                    elif update and not cmap_input:
                        # set alpha value for transparency
                        cmap = plt.cm.Reds
                        my_cmap = cmap(np.arange(cmap.N))
                        # Set alpha
                        my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
                        # Create new colormap
                        my_cmap = ListedColormap(my_cmap)
                        self.colormap = my_cmap
                    else:
                        # set alpha value for transparency
                        cmap = plt.get_cmap(cmap_input)
                        my_cmap = cmap(np.arange(cmap.N))
                        # Set alpha
                        my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
                        # Create new colormap
                        my_cmap = ListedColormap(my_cmap)
                        self.colormap = my_cmap
                else:
                    if not update:
                        self.colormap = cm.PuBu;
                    elif update and not cmap_input:
                        # set alpha value for transparency
                        cmap = plt.cm.PuBu
                        my_cmap = cmap(np.arange(cmap.N))
                        # Set alpha
                        my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
                        # Create new colormap
                        my_cmap = ListedColormap(my_cmap)
                        self.colormap = my_cmap
                    else:
                        # set alpha value for transparency
                        cmap = plt.get_cmap(cmap_input)
                        my_cmap = cmap(np.arange(cmap.N))
                        # Set alpha
                        my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
                        # Create new colormap
                        my_cmap = ListedColormap(my_cmap)
                        self.colormap = my_cmap
                clim[0] = np.amin(self.slice)
                clim[1] = np.amax(self.slice)
                if (clim[0] == 0.0) and (clim[1] == 0.0):
                    clim[0] = 0.0
                    clim[1] = 1.0
            else:
                if not update:
                    self.colormap = cm.RdGy;
                elif update and not cmap_input:
                    # set alpha value for transparency
                    cmap = plt.cm.RdGy
                    my_cmap = cmap(np.arange(cmap.N))
                    # Set alpha
                    my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
                    # Create new colormap
                    my_cmap = ListedColormap(my_cmap)
                    self.colormap = my_cmap
                else:
                    # set alpha value for transparency
                    cmap = plt.get_cmap(cmap_input)
                    my_cmap = cmap(np.arange(cmap.N))
                    # Set alpha
                    my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
                    # Create new colormap
                    my_cmap = ListedColormap(my_cmap)
                    self.colormap = my_cmap
                clim[0] = -np.amax(abs(self.slice))
                clim[1] = np.amax(abs(self.slice))              
                if (clim[0] == 0.0) and (clim[1] == 0.0):
                    clim[0] = -1.0
                    clim[1] = 1.0

        elif self.g3d.get_type() == pp_defs.hipace.h5.g3dtypes.field:
            if not update:
                self.colormap = cm.RdBu;
            elif update and not cmap_input:
                # set alpha value for transparency
                cmap = plt.cm.RdBu
                my_cmap = cmap(np.arange(cmap.N))
                # Set alpha
                alpha_array = np.linspace(1, 0, cmap.N/2)
                alpha_array = np.append(alpha_array, np.linspace(0, 1, cmap.N/2) )
                my_cmap[:,-1] = alpha_array #np.linspace(0, 1, cmap.N)
                # Create new colormap
                my_cmap = ListedColormap(my_cmap)
                self.colormap = my_cmap
            else:
                # set alpha value for transparency
                cmap = plt.get_cmap(cmap_input)
                my_cmap = cmap(np.arange(cmap.N))
                # Set alpha
                my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
                # Create new colormap
                my_cmap = ListedColormap(my_cmap)
                self.colormap = my_cmap
            clim[0] = -np.amax(abs(self.slice))
            clim[1] = np.amax(abs(self.slice))
            if (clim[0] == 0.0) and (clim[1] == 0.0):
                clim[0] = -1.0
                clim[1] = 1.0            
        elif self.g3d.get_type() == pp_defs.hipace.h5.g3dtypes.current:
            if not update:
                self.colormap = cm.PuOr;
            elif update and not cmap_input:
                # set alpha value for transparency
                cmap = plt.cm.PuOr
                my_cmap = cmap(np.arange(cmap.N))
                # Set alpha
                my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
                # Create new colormap
                my_cmap = ListedColormap(my_cmap)
                self.colormap = my_cmap
            else:
                # set alpha value for transparency
                cmap = plt.get_cmap(cmap_input)
                my_cmap = cmap(np.arange(cmap.N))
                # Set alpha
                my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
                # Create new colormap
                my_cmap = ListedColormap(my_cmap)
                self.colormap = my_cmap
            clim[0] = -np.amax(abs(self.slice))
            clim[1] = np.amax(abs(self.slice))
            if (clim[0] == 0.0) and (clim[1] == 0.0):
                clim[0] = -1.0
                clim[1] = 1.0
        else:
            print('no hipace type data detected!')
                
        if self.args.diff:
            self.colormap = cm.RdBu

        if self.args.clog:
            self.app_str += '_log'
            min_nonzero = np.min(abs(self.slice[np.where(self.slice != 0.0)]))
            self.slice[np.where(self.slice == 0)] = min_nonzero
            self.slice = abs(self.slice)
            self.colormap = cm.Greys
            clim[0] = min_nonzero
            clim[1] = np.max(self.slice)

        if self.args.clim != None and not update:
            clim = list(self.args.clim)
        if clim_input != None:
            clim = list(clim_input)

        self.clim = clim

    def plot( self, ifsave=True ):
        if self.args.verbose: 
            sys.stdout.write('Generating slice plot\n')
            sys.stdout.flush()
        saveformat = self.args.file_format
        filesuffix = '_%06.f' % (np.floor(self.g3d.get_time()))

        if self.args.save_prefix != parsedefs.save.prefix:
            fileprefix = self.args.save_prefix
        else:
            fileprefix = self.g3d.get_name()

        if self.g3d.is_subgrid():
            self.app_str += '_subgrid'


        fig = plt.figure()
        if self.args.ptype == 'pcolormesh':
            cax = plt.pcolormesh(self.x_array,
                                 self.y_array,
                                 self.slice,
                                 vmin=self.clim[0], vmax=self.clim[1])
        elif self.args.ptype == 'contourf':
            levels = MaxNLocator(nbins='512', steps=[1, 2, 4, 5, 10]).tick_values(self.clim[0], self.clim[1])
            
            #selecting correct extend method
            if np.amin(self.slice) < self.clim[0] and np.amax(self.slice) > self.clim[1]:
                extend = 'both'
            elif np.amin(self.slice) < self.clim[0] and np.amax(self.slice) <= self.clim[1]:
                extend = 'min'
            elif np.amin(self.slice) >= self.clim[0] and np.amax(self.slice) > self.clim[1]:
                extend = 'max'
            elif np.amin(self.slice) >= self.clim[0] and np.amax(self.slice) <= self.clim[1]:
                extend = 'neither'
            else:
                print('Error: unexpected case, couldn\'t extend in the correct way!')
                extend = 'neither'
            cax = plt.contourf( self.x_array,
                                self.y_array,
                                self.slice,
                                levels=levels,
                                vmin=self.clim[0], vmax=self.clim[1],
                                cmap=self.colormap,
                                extend=extend)
        elif self.args.ptype == 'pcolor':
            cax = plt.pcolor( self.x_array,
                                self.y_array,
                                self.slice,
                                vmin=self.clim[0], vmax=self.clim[1])            
        elif self.args.ptype == 'imshow':
            sys.stdout.write('ERROR: imshow not implemented yet!\n')
            sys.stdout.flush()
            sys.exit(1)  
        elif self.args.ptype == 'pcolorfast':    
            sys.stdout.write('ERROR: pcolorfast not implemented yet!\n')
            sys.stdout.flush()
            sys.exit(1)     

        cax.cmap = self.colormap

        if self.args.clog:
            cax.norm = matplotlib.colors.LogNorm(vmin=self.clim[0], vmax=self.clim[1])

        ax = plt.gca()
        ax.tick_params(labelsize=self.args.axticklabelsize)
        ax.set_ylabel(self.ylabel, fontsize=self.args.fontsize)
        ax.set_xlabel(self.xlabel, fontsize=self.args.fontsize)
        cbar = fig.colorbar(cax)
        cbar.ax.set_ylabel( gen_pretty_grid_name( self.g3d.get_name() ), fontsize=self.args.fontsize )
        
        #manually setting cbar ticks to avoid cutoff of the last tick
        ticks = MaxNLocator().tick_values(self.clim[0], self.clim[1])
        cbar.set_ticks ( ticks )
        
        if self.args.xlim != None:
            plt.xlim(self.args.xlim[0], self.args.xlim[1])
        if self.args.ylim != None:
            plt.ylim(self.args.ylim[0], self.args.ylim[1])

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

        if self.args.savepath == None:
            savepath = self.root_savepath + self.__relsavepath
        else:
            savepath = self.args.savepath + self.__relsavepath

        mkdirs_if_nexist(savepath)


        savename = "%s_%s%s%s.%s" % (  fileprefix, \
                                       self.args.plane, \
                                       self.app_str, \
                                       filesuffix, \
                                       saveformat )

        if saveformat==parsedefs.file_format.png:
            fig.savefig( savepath + '/' + savename,
                      format=saveformat,
                      dpi=300)
        else:
            fig.savefig(  savepath + '/' + savename,
                          format=saveformat)
        if self.args.verbose: 
            sys.stdout.write('Saved "%s" at: %s\n' % (savename,savepath))
            sys.stdout.flush()

        if self.args.ifshow: plt.show()
        plt.close(fig)
        
    def update_fig( self, ifsave=True, update=False,  fig=None, diff=False, diff_data=None, counts=0, line_data=False ):
        if self.args.verbose: 
            sys.stdout.write('Generating slice plot\n')
            sys.stdout.flush()

        if self.g3d.is_subgrid():
            self.app_str += '_subgrid'
        
        if diff:
            self.slice = self.slice - diff_data
        
        #selecting correct extend method
        if np.amin(self.slice) < self.clim[0] and np.amax(self.slice) > self.clim[1]:
            extend = 'both'
        elif np.amin(self.slice) < self.clim[0] and np.amax(self.slice) <= self.clim[1]:
            extend = 'min'
        elif np.amin(self.slice) >= self.clim[0] and np.amax(self.slice) > self.clim[1]:
            extend = 'max'
        elif np.amin(self.slice) >= self.clim[0] and np.amax(self.slice) <= self.clim[1]:
            extend = 'neither'
        else:
            print('Error: unexpected case, couldn\'t extend in the correct way!')
            extend = 'neither'
        
        if not update:
            fig = plt.figure()
        if self.args.ptype == 'pcolormesh':
            cax = plt.pcolormesh(self.x_array,
                                 self.y_array,
                                 self.slice,
                                 vmin=self.clim[0], vmax=self.clim[1],
                                 cmap=self.colormap)
        elif self.args.ptype == 'contourf':
            levels = MaxNLocator(nbins='512', steps=[1, 2, 4, 5, 10]).tick_values(self.clim[0], self.clim[1])
            

            cax = plt.contourf( self.x_array,
                                self.y_array,
                                self.slice,
                                levels=levels,
                                vmin=self.clim[0], vmax=self.clim[1],
                                cmap=self.colormap,
                                extend=extend)
        elif self.args.ptype == 'pcolor':
            cax = plt.pcolor( self.x_array,
                                self.y_array,
                                self.slice,
                                cmap=self.colormap,
                                vmin=self.clim[0], vmax=self.clim[1])            
        elif self.args.ptype == 'imshow':
            sys.stdout.write('ERROR: imshow not implemented yet!\n')
            sys.stdout.flush()
            sys.exit(1)  
        elif self.args.ptype == 'pcolorfast':    
            sys.stdout.write('ERROR: pcolorfast not implemented yet!\n')
            sys.stdout.flush()
            sys.exit(1)     

        #cax.cmap = self.colormap
        cbarmap = plt.cm.ScalarMappable(cmap=self.colormap)
        cbarmap.set_array(self.clim[1]-self.clim[0])

        # Note on colorbar: boundaries have to be set manually, because otherwise there will be ugly stripes
        # afterwards the ticks have to set manually as well, set them at the correct place
        cbarmap.set_clim(self.clim[0], self.clim[1])
        
        print(line_data)
        if line_data:
            offset = 0.15
        else:
            offset = 0

        if self.g3d.get_type() == pp_defs.hipace.h5.g3dtypes.field:
            cbar = plt.colorbar(cbarmap, boundaries=np.arange(self.clim[0]-0.0002,self.clim[1]+0.0002,0.0001), extend=extend, pad=0.02+offset, fraction=0.08+4*counts**(1/2)/100.0) #, pad=0.1 )
        else: 
            cbar = plt.colorbar(cbarmap, boundaries=np.arange(self.clim[0],self.clim[1]+0.0002,0.0001), extend=extend, pad=0.02+offset, fraction=0.08+4*counts**(1/2)/100.0) #, pad=0.1 )
        ticks = MaxNLocator().tick_values(self.clim[0], self.clim[1])
        cbar.set_ticks ( ticks )
        cbar.set_clim([self.clim[0], self.clim[1]])
        
        if self.args.clog:
            cax.norm = matplotlib.colors.LogNorm(vmin=self.clim[0], vmax=self.clim[1])

        ax = plt.gca()
        ax.tick_params(labelsize=self.args.axticklabelsize)
        ax.set_ylabel(self.ylabel, fontsize=self.args.fontsize, labelpad=self.args.ylabelpad)
        ax.set_xlabel(self.xlabel, fontsize=self.args.fontsize)
        #cbar = fig.colorbar(cax)
        cbar.ax.set_title( gen_pretty_grid_name( self.g3d.get_name() ), fontsize=self.args.fontsize, pad=self.args.cbarpad )
        cbar.ax.tick_params(labelsize=self.args.axticklabelsize)
        #manually setting cbar ticks to avoid cutoff of the last tick
        ticks = MaxNLocator(9).tick_values(self.clim[0], self.clim[1])
        cbar.set_ticks ( ticks )
        
        if self.args.xlim != None:
            plt.xlim(self.args.xlim[0], self.args.xlim[1])
        if self.args.ylim != None:
            plt.ylim(self.args.ylim[0], self.args.ylim[1])

        # axpos1 = ax.get_position() # get the original position 
        # cbaxpos1 = cbar.ax.get_position() # get the original position 
        # xshift = 0.05;
        # wcompress = 1.05
        # axpos2 = [axpos1.x0 + xshift, axpos1.y0,  (axpos1.width - xshift) / wcompress, axpos1.height ]
        # ax.set_position(axpos2) # set a new position 
        # axpos2 = ax.get_position()
        # cbaxpos2 = [cbaxpos1.x0 + (axpos2.x0 - axpos1.x0) + (axpos2.width - axpos1.width), \
        #             cbaxpos1.y0, \
        #             cbaxpos1.width, \
        #             cbaxpos1.height ] 
        # 
        # cbar.ax.set_position(cbaxpos2)
        # 
        return fig
        
    def save_fig(self, fig):
        
        saveformat = self.args.file_format
        filesuffix = '_%06.f' % (np.floor(self.g3d.get_time()))

        if self.args.save_prefix != parsedefs.save.prefix:
            fileprefix = self.args.save_prefix
        else:
            fileprefix = self.g3d.get_name()
        
        if self.args.savepath == None:
            savepath = self.root_savepath + self.__relsavepath
        else:
            savepath = self.args.savepath + self.__relsavepath

        mkdirs_if_nexist(savepath)


        savename = "%s_%s%s%s.%s" % (  fileprefix, \
                                       self.args.plane, \
                                       self.app_str, \
                                       filesuffix, \
                                       saveformat )

        if saveformat==parsedefs.file_format.png:
            fig.savefig( savepath + '/' + savename,
                      format=saveformat,
                      dpi=300)
        else:
            fig.savefig(  savepath + '/' + savename,
                          format=saveformat)
        if self.args.verbose: 
            sys.stdout.write('Saved "%s" at: %s\n' % (savename,savepath))
            sys.stdout.flush()

        if self.args.ifshow: plt.show()
        #fig.clear()
        plt.close(fig)


class G3d_plot_line(G3d_plot):
    def __init__(self, file, args):
        G3d_plot.__init__(self, file, args)

        self.__relsavepath = '/g3d-line'

        self.set_xaxis( self.args.lineax )
        if self.args.verbose: 
            sys.stdout.write('Reading data...\n')
            sys.stdout.flush()
        self.read()
        if self.args.verbose: 
            sys.stdout.write('Read-in completed.\n')
            sys.stdout.flush()
        self.set_yaxis()


    def set_yaxis( self ):
        if self.args.if_integrate:
            self.ylabel = r'$k_p^2 \int \int$' + gen_pretty_grid_name( self.g3d.get_name() )
            if self.args.lineax == 'z':
                self.ylabel = self.ylabel + r'$\,dx dy$'
            if self.args.lineax == 'x':
                self.ylabel = self.ylabel + r'$\,dy dz$'
            if self.args.lineax == 'y':
                self.ylabel = self.ylabel + r'$\,dx dz$'
        else:
            self.ylabel = gen_pretty_grid_name( self.g3d.get_name() )
        ylim = [0.0, 0.0]
        # define axis labels and arrays
        if self.g3d.get_type() == pp_defs.hipace.h5.g3dtypes.density:
            ylim[0] = np.amin(self.line)
            ylim[1] = np.amax(self.line)
        elif self.g3d.get_type() == pp_defs.hipace.h5.g3dtypes.field:
            self.colormap = cm.coolwarm
            ylim[0] = -np.amax(abs(self.line))
            ylim[1] = np.amax(abs(self.line))
        elif self.g3d.get_type() == pp_defs.hipace.h5.g3dtypes.current:
            self.colormap = cm.coolwarm
            ylim[0] = -np.amax(abs(self.line))
            ylim[1] = np.amax(abs(self.line))
        self.ylim = ylim

    def read( self ):
        # read line
        if self.args.lout_idx != None:
            lout_idx = list(self.args.lout_idx)

        if 'z' == self.args.lineax:
            if self.args.if_integrate:
                self.line = self.g3d.read_integrate(ax1=True,ax2=True)
            elif self.args.avgax != None:
                if self.args.avgax == 'x':
                    self.line = self.g3d.read_avgx(dim=1,ax2=True)
                elif self.args.avgax == 'y':
                    self.line = self.g3d.read_avgx(dim=2,ax1=True)
                else:
                    sys.stdout.write('ERROR: Cannot average along plotted axis!\n')
                    sys.stdout.flush()
                    sys.exit(1)
            else:
                if self.args.lout_idx == None:
                    # Default: central lineout
                    if self.args.piccode == pp_defs.code.inferno:
                        self.line= self.g3d.read(i1 = self.g3d.get_nx(1)/2)
                        print(np.shape(self.line))
                    else:
                        idx1 = math.floor(self.g3d.get_nx(1)/2) - 1
                        idx2 = math.floor(self.g3d.get_nx(2)/2) - 1
                        self.line = self.g3d.read(i1=idx1, i2=idx2)
                        if self.g3d.get_nx(1)%2 == 0 and self.g3d.get_nx(2)%2 == 0:
                            line01 = self.g3d.read(i1=idx1, i2=idx2+1)
                            line10 = self.g3d.read(i1=idx1+1, i2=idx2)
                            line11 = self.g3d.read(i1=idx1+1, i2=idx2+1)
                            self.line = ( self.line + line01 + line10 + line11 )/4
                        elif self.g3d.get_nx(1)%2 == 1 and self.g3d.get_nx(2)%2 == 0:
                            self.line = ( self.line + self.self.g3d.read(i1=idx1, i2=idx2+1) )/2
                        elif self.g3d.get_nx(1)%2 == 0 and self.g3d.get_nx(2)%2 == 1:
                            self.line = ( self.line + self.g3d.read(i1=idx1+1, i2=idx2) )/2
                else:
                    if self.args.piccode == pp_defs.code.inferno:
                        self.line= self.g3d.read(i1 = lout_idx[0])
                    else:
                        self.line = self.g3d.read(i1=lout_idx[0], i2=lout_idx[1])

        elif 'x' == self.args.lineax:
            if self.args.if_integrate:
                self.line = self.g3d.read_integrate(ax0=True,ax2=True)
            elif self.args.avgax != None:
                if self.args.avgax == 'y':
                    self.line = self.g3d.read_avgx(dim=2,ax0=True)
                elif self.args.avgax == 'z':
                    self.line = self.g3d.read_avgx(dim=0,ax2=True)
                else:
                    sys.stdout.write('ERROR: Cannot average along plotted axis!\n')
                    sys.stdout.flush()
                    sys.exit(1)                 
            else:
                if (self.args.lout_idx == None) and (self.args.lout_zeta_pos == None):
                    if self.args.piccode == pp_defs.code.inferno:
                        self.line= self.g3d.read(i0 = g3d.get_nx(1)/2)
                    else:
                        # Default: central lineout
                        idx1 = math.floor(g3d.get_nx(0)/2) - 1
                        idx2 = math.floor(g3d.get_nx(2)/2) - 1
                        self.line = self.g3d.read(i0=idx1, i2=idx2)
                        if g3d.get_nx(0)%2 == 0 and g3d.get_nx(2)%2 == 0:
                            line01 = self.g3d.read(i0=idx1, i2=idx2+1)
                            line10 = self.g3d.read(i0=idx1+1, i2=idx2)
                            line11 = self.g3d.read(i0=idx1+1, i2=idx2+1)
                            self.line = ( self.line + line01 + line10 + line11 )/4
                        elif g3d.get_nx(0)%2 == 1 and g3d.get_nx(2)%2 == 0:
                            self.line = ( self.line + self.g3d.read(i0=idx1, i2=idx2+1) )/2
                        elif g3d.get_nx(0)%2 == 0 and g3d.get_nx(2)%2 == 1:
                            self.line = ( self.line + self.g3d.read(i0=idx1+1, i2=idx2) )/2
                elif (self.args.lout_idx != None) and (self.args.lout_zeta_pos == None):
                    self.line = self.g3d.read(i0=lout_idx[0], i2=lout_idx[1])
                elif (self.args.lout_idx == None) and (self.args.lout_zeta_pos != None):
                    if self.args.piccode == pp_defs.code.inferno:
                        self.line= self.g3d.read(x0=self.args.lout_zeta_pos)
                    else:
                        self.line = self.g3d.read(x0=self.args.lout_zeta_pos, x2=0.0)
                else:
                    sys.stdout.write('ERROR: lineout-index can''t be used in conjunction with lineout-zeta-position!\n')
                    sys.stdout.flush()
                    sys.exit(1)                

        elif 'y' == self.args.lineax:
            if self.args.if_integrate:
                self.line = self.g3d.read_integrate(ax0=True,ax1=True)
            elif self.args.avgax != None:
                if self.args.avgax == 'x':
                    self.line = self.g3d.read_avgx(dim=1,ax0=True)
                elif self.args.avgax == 'z':
                    self.line = self.g3d.read_avgx(dim=0,ax1=True)
                else:
                    sys.stdout.write('ERROR: Cannot average along plotted axis!\n')
                    sys.stdout.flush()
                    sys.exit(1)                 
            else:                
                if (self.args.lout_idx == None) and (self.args.lout_zeta_pos == None):
                    # Default: central lineout
                    idx1 = math.floor(self.g3d.get_nx(0)/2) - 1
                    idx2 = math.floor(self.g3d.get_nx(1)/2) - 1
                    self.line = self.g3d.read(i0=idx1, i1=idx2)
                    if self.g3d.get_nx(0)%2 == 0 and self.g3d.get_nx(1)%2 == 0:
                        line01 = self.g3d.read(i0=idx1, i1=idx2+1)
                        line10 = self.g3d.read(i0=idx1+1, i1=idx2)
                        line11 = self.g3d.read(i0=idx1+1, i1=idx2+1)
                        self.line = ( self.line + line01 + line10 + line11 )/4
                    elif self.g3d.get_nx(0)%2 == 1 and self.g3d.get_nx(1)%2 == 0:
                        self.line = ( self.line + self.g3d.read(i0=idx1, i1=idx2+1) )/2
                    elif self.g3d.get_nx(0)%2 == 0 and self.g3d.get_nx(1)%2 == 1:
                        self.line = ( self.line + self.g3d.read(i0=idx1+1, i1=idx2) )/2
                elif (self.args.lout_idx != None) and (self.args.lout_zeta_pos == None):
                    self.line = self.g3d.read(i0=lout_idx[0], i1=lout_idx[1])
                elif (self.args.lout_idx == None) and (self.args.lout_zeta_pos != None):
                    self.line = self.g3d.read(x0=self.args.lout_zeta_pos, x1=0.0)
                else:
                    sys.stdout.write('ERROR: lineout-index can''t be used in conjunction with lineout-zeta-position!\n')
                    sys.stdout.flush()
                    sys.exit(1)                                               

        if self.is_number_density():
            self.line = np.abs(self.line)

    def plot( self, ifsave=True, fig=False, smooth=False ):
        if self.args.verbose: 
            sys.stdout.write('Generating line plot...\n')
            sys.stdout.flush()
        saveformat = self.args.file_format
        filesuffix = '_%06.f' % (np.floor(self.g3d.get_time()))

        if self.args.save_prefix != parsedefs.save.prefix:
            fileprefix = self.args.save_prefix
        else:
            fileprefix = self.g3d.get_name()

        if self.g3d.is_subgrid():
            self.app_str += '_subgrid_' 

        if fig:
            self.app_str += '_overlay'

        if self.args.if_integrate:
            self.app_str += '_int'
        elif self.args.avgax != None:
            self.app_str += '_avg' + self.args.avgax

        if not fig:
            fig = plt.figure()
            
        ax = plt.gca()
        if fig:
            ax = ax.twinx()
            ax.set_xlim(min(self.x_array), max(self.x_array))
        else:
            ax.set_xlabel(self.xlabel, fontsize=self.args.fontsize)
        ax.tick_params(labelsize=self.args.axticklabelsize)    
        ax.set_ylabel(self.ylabel, fontsize=self.args.fontsize, labelpad=self.args.y2labelpad)
        
        
        if self.args.absylog:
            self.app_str += '_absylog'
            nonzero_idx = np.where( abs(self.line) > 0.0 )
            plt.semilogy( self.x_array[nonzero_idx],
                          abs(self.line[nonzero_idx]))
        elif smooth:
            plt.plot( self.x_array,
                      gaussian_filter1d(self.line, smooth))
        else:
            plt.plot( self.x_array,
                      self.line)
            

        

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

        if self.args.xlim != None:
            plt.xlim(self.args.xlim[0], self.args.xlim[1])
        if self.args.ylim != None:
            plt.ylim(self.args.ylim[0], self.args.ylim[1])

        if self.args.ylim2:
            plt.ylim(self.args.ylim2[0], self.args.ylim2[1])
            ax.set_ylim( self.args.ylim2[0], self.args.ylim2[1] )

        savename = "%s_%s%s%s" % (  fileprefix, \
                                   self.args.lineax, \
                                   self.app_str, \
                                   filesuffix)

        #plt.title(savename)

        if self.args.savepath == None:
            savepath = self.root_savepath + self.__relsavepath
        else:
            savepath = self.args.savepath + self.__relsavepath

        mkdirs_if_nexist(savepath)
        if self.args.smoothing:
            savename += '_smoothed'
            
        if saveformat==parsedefs.file_format.png:
            fig.savefig( savepath + '/' + savename + '.' + saveformat,
                      format=saveformat,
                      dpi=self.args.dpi)
        else:
            fig.savefig( savepath + '/' + savename + '.' + saveformat,
                          format=saveformat)
        if self.args.verbose:
            sys.stdout.write('Saved "%s.%s" at: %s\n' % (savename, saveformat, savepath ))
            sys.stdout.flush()

        if not self.args.h5plot_off: 
            h5lp = H5Plot()
            h5lp.inherit_matplotlib_line_plots(ax)
            h5lp.write(savepath + '/' + savename + '.h5')

        if self.args.ifshow: plt.show()
        plt.close(fig)


def slice(args):

    h5flist = H5FList(args.path, h5ftype='g3d')
    flist = h5flist.get()

    if args.line_data:
        line_data = True
    else:
        line_data = False
    for file in flist:
        timestamp = file.split("_")[-1].split(".h5")[0]
        
        g3d_p = G3d_plot_slice(file, args)
        if not args.data2 and not args.data3 and not args.line_data:
            g3d_p.plot()
        elif args.diff:
            for data2_files in args.data2:
                if timestamp in data2_files or args.manual:
                    g3d_p2 = G3d_plot_slice(data2_files, args)
                    print(np.shape(g3d_p2.slice))
                    fig = g3d_p.update_fig(diff=True, diff_data=g3d_p2.slice, line_data=args.line_data)
                    g3d_p.save_fig(fig)
                    plt.close(fig)
        else:
            if not args.diff:
                fig = g3d_p.update_fig(line_data=args.line_data)
                if args.data2:
                    fig.set_size_inches(8.3, 6)
                    for data2_files in args.data2:
                        if timestamp in data2_files or args.manual:
                            g3d_p2 = G3d_plot_slice(data2_files, args)
                            g3d_p2.set_cmap( update=True, clim_input=args.clim2, cmap_input=args.cm2 )
                            fig = g3d_p2.update_fig( fig=fig, update=True, counts=1, line_data=args.line_data)
                            if args.data3:
                                fig.set_size_inches(8.6, 6)
                                for data3_files in args.data3:
                                    if timestamp in data3_files or args.manual:
                                        g3d_p3 = G3d_plot_slice(data3_files, args)
                                        g3d_p3.set_cmap( update=True, clim_input=args.clim3, cmap_input=args.cm3 )
                                        fig = g3d_p3.update_fig( fig=fig, update=True, counts=2, line_data=args.line_data)
                if args.line_data:
                    for line_data_files in args.line_data:
                        if timestamp in line_data_files or args.manual:
                            g3d_pline = G3d_plot_line(line_data_files, args)
                            g3d_pline.set_yaxis()
                            g3d_pline.plot(fig=fig, smooth=args.smoothing)
                else:
                    g3d_p.save_fig(fig)
                    plt.close(fig)


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

    if args.latexfont:
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif') 

    args.func(args)


if __name__ == "__main__":
    main()         