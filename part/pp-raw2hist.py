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
from matplotlib import cm

mypath = os.path.dirname(os.path.abspath( __file__ ))
incpath = os.path.split(mypath)[0] + '/inc'
sys.path.append(incpath)
import pp_defs
from pp_h5dat import HiRAW
from pp_h5dat import H5FList
from pp_plt_tools import saveas_png
from pp_plt_tools import saveas_eps_pdf
from pp_raw_ana import TagSelect


def raw2hist_parser():

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
                          help = "Print info (Default).")
    parser.add_argument(  "-q", "--quiet",
                          dest = "verbose",
                          action = "store_false",
                          help = "Don't print info.")
    parser.add_argument(  "--show",
                          dest = "ifshow",
                          action = "store_true",
                          default = False,
                          help = "Show figure.")
    parser.add_argument(  '--zeta-range',
                          help='zeta range',
                          action='store',
                          dest="zeta_range",
                          metavar=('ZETA_MIN', 'ZETA_MAX'),
                          nargs=2,
                          type=float,
                          default=None)
    parser.add_argument(  "--latexfont",
                          dest = "latexfont",
                          action="store_true",
                          default=False,
                          help = "Use LaTeX font (Default: %(default)s).")  
    parser.add_argument(  "--psel",
                          dest = "part_selection_criterion",
                          action="store",
                          default='',
                          type=str,
                          help = "Particle selection criterion (e.g. 'x**2 + y**2 < 0.1').")
    parser.add_argument(  "--tfout",
                          dest = "tagfile_out",
                          action="store",
                          default=None,
                          type=str,
                          help = "Save selected particle tags to specified tag file (Default: %(default)s).")    
    parser.add_argument(  "--tfin",
                          dest = "tagfile_in",
                          action="store",
                          default=None,
                          type=str,
                          help = "Particle selection from specified tag file (Default: %(default)s).")
    parser.add_argument(  "--stride",
                          dest = "stride",
                          action="store",
                          default=1,
                          type=int,
                          help = "For a selected filelist, process every 'stride' file (Default: %(default)s).")                                                                                                          
    return parser


def r2h_1d_subparser(subparsers, parent_parser):
    parser = subparsers.add_parser( "1d", parents=[parent_parser],
                                    help="1D histogram")    
    # Line plot specific arguments
    parser.add_argument(  '--psv', '--phase-space-var',
                          action='store',
                          dest="psv",
                          metavar="VAR",
                          choices=[ 'x', 'y', 'z', 'r', 'px', 'py', 'pz'],
                          default=None,
                          required=True,
                          help= 'Phase space variable')
    parser.add_argument(  "-s", "--save-path",
                          action="store",
                          dest="savepath",
                          metavar="PATH",
                          default='plots/hist1d',
                          help = """Path to which generated files will be saved.
                              (Default: './')""")
    parser.add_argument(  "-f", "--format",
                          action='store',
                          dest="file_format",
                          metavar="FORMAT",
                          choices=[ 'png',
                                    'pdf',
                                    'eps',],
                          default='eps',
                          help= """Format of output file (Default: eps).""")
    parser.add_argument(  '-r','--range',
                          action='store',
                          dest="range",
                          metavar=('XMIN', 'XMAX'),
                          nargs=2,
                          type=float,
                          default=None,
                          help='Range of histogram.')
    parser.add_argument(  '-n','--nbins',
                          help='Number of bins.',
                          action='store',
                          dest="nbins",
                          metavar='NBINS',
                          type=int,
                          default=None)                              
    return parser

def r2h_2d_subparser(subparsers, parent_parser):
    parser = subparsers.add_parser( "2d", parents=[parent_parser],
                                    help="2d histogram")    
    # Slice plot specific arguments
    parser.add_argument(  '--psv', '--phase-space-vars',
                          help='Phase space variables',
                          action='store',
                          dest="psv",
                          metavar=('VAR1', 'VAR2'),
                          nargs=2,
                          required=True,
                          default=None)
    parser.add_argument(  "--cscale",
                          action='store',
                          dest="cscale",
                          metavar="CSCALE",
                          choices=[ "lin", "log",],
                          default="lin",
                          help= "z-axis type (Default: %(default)s).")
    parser.add_argument(  '--cblim',
                          help='Colorbar axis limits',
                          action='store',
                          dest="cblim",
                          metavar=('CBMIN', 'CBMAX'),
                          nargs=2,
                          type=float,
                          default=None)
    parser.add_argument(  '-n', '--nbins',
                          help='Number of bins in each dimension.',
                          action='store',
                          dest="nbins",
                          metavar=('NBINS1', 'NBINS2'),
                          nargs=2,
                          type=int,
                          default=None)
    parser.add_argument(  '--range',
                          help='Range in both dimensions.',
                          action='store',
                          dest="range",
                          metavar=('XMIN', 'XMAX', 'YMIN', 'YMAX'),
                          nargs=4,
                          type=float,
                          default=None)                                                        
    parser.add_argument(  "-s", "--save-path",
                          action="store",
                          dest="savepath",
                          metavar="PATH",
                          default='plots/hist2d',
                          help = """Path to which generated files will be saved.
                              (Default: './')""")
    parser.add_argument(  "-f", "--format",
                          action='store',
                          dest="file_format",
                          metavar="FORMAT",
                          choices=[ 'png',
                                    'pdf',
                                    'eps',],
                          default='png',
                          help= """Format of output file (Default: png).""")    
    return parser


def get_props(raw,psv_str):
    if psv_str == 'x':
       label = r'$k_p x$' 
       psv = raw.get('x2')
       savename = 'x'
    elif psv_str == 'y':
       label = r'$k_p y$' 
       psv = raw.get('x3')
       savename = 'y'       
    elif psv_str == 'z' or psv_str == 'zeta':
       label = r'$k_p \zeta$' 
       psv = raw.get('x1')
       savename = 'zeta'
    elif psv_str == 'r':
       label = r'$k_p r$' 
       psv = np.sqrt(np.power(raw.get('x2'),2) + np.power(raw.get('x3'),2)) 
       savename = 'r'    
    elif psv_str == 'px':
       label = r'$p_x/m c$' 
       psv = raw.get('p2')
       savename = 'px'
    elif psv_str == 'py':
       label = r'$p_y/m c$' 
       psv = raw.get('p3')
       savename = 'py'
    elif psv_str == 'pz':
       label = r'$p_z/m c$' 
       psv = raw.get('p1')
       savename = 'pz'
    else:
        print("ERROR: Phase space variable '%s' not recognized!" % psv_str)
        sys.exit(1)
    return psv, label, savename


def get_zeta_range_str(zeta_range):
    if zeta_range != None:
        zr_str = ('_zeta_range_%0.2f_%0.2f' % (zeta_range[0], zeta_range[1]))
    else:
        zr_str = ''
    return zr_str 


def oneD(raw, args):

    appstr = ''

    psv, xlabel, savename = get_props(raw,args.psv)

    if args.nbins == None:
        nbins = np.int( 10 + np.power(np.log(raw.get_npart()), 3) )     
    else:
        nbins = args.nbins

    hist, bin_edges = np.histogram( psv,\
                                    bins=nbins,\
                                    weights=raw.get('q'),\
                                    density=True,\
                                    range=args.range)
    x_array = bin_edges[0:-1] + (bin_edges[1] - bin_edges[0])/2

    if args.psv == 'r':
        ## Scaling dr with 1/r
        # rmax = np.max(psv)
        # r_range = [0, rmax]
        # r_lsp = np.linspace(0,rmax,nbins)
        # lsp_centers = r_lsp[:-1] + np.diff(r_lsp)
        # dr_unnormed = np.divide(np.diff(r_lsp),lsp_centers)
        # dr_normed =  rmax * dr_unnormed/np.sum(dr_unnormed) 
        # r_bin_edges = np.insert(np.cumsum(dr_normed),0,0)
        #hist, bin_edges = np.histogram(psv,bins=r_bin_edges,weights=raw.get('q'),density=True,range=args.range)

        ## Scaling number of hits with 1/r
        hist = np.divide(hist,x_array)/(2*math.pi)

    fig = plt.figure(figsize=(6,5))
    plt.plot( x_array, hist )
    ax = plt.gca()    
    ax.set_xlabel(xlabel, fontsize=14)
    plt.gcf().subplots_adjust(left=0.15, bottom=0.15)   

    if args.part_selection_criterion != '' or args.tagfile_in != None:
        appstr += '_sel'
    
    raw.read_attrs()
    timestamp = '_%08.1f' % (raw.get_time()) 
    savename += get_zeta_range_str(args.zeta_range)  \
                + appstr \
                + timestamp
    
    if args.file_format == 'png':
        saveas_png(fig, args.savepath, savename, verbose=True)
    else:
        saveas_eps_pdf(fig, args.savepath, savename, h5plot=True, verbose=True, fformat=args.file_format)
    plt.close(fig)


def twoD(raw, args):

    appstr = ''

    varx, xlabel, savenamex = get_props(raw,args.psv[0])
    vary, ylabel, savenamey = get_props(raw,args.psv[1])

    if args.nbins == None:
        nbins = np.int( 10 + np.power(np.log(raw.get_npart()), 2.2) )
    else:
        nbins = args.nbins

    if args.range != None:
        hrange = [[args.range[0], args.range[1]], [args.range[2], args.range[3]]]
    else:
        hrange = None

    if args.verbose:
        print("Generating 2d histogram...")
    H, xedges, yedges = np.histogram2d(varx,vary,bins=nbins,weights=raw.get('q'),range=hrange)
    x_array = xedges[0:-1] + (xedges[1] - xedges[0])/2
    y_array = yedges[0:-1] + (yedges[1] - yedges[0])/2

    if args.verbose:
        print("Generating plot...")
    fig = plt.figure(figsize=(6,5))
    cax = plt.pcolor(x_array,
                     y_array,
                     H.transpose(), cmap='PuBu') 
    plt.gcf().subplots_adjust(left=0.15, bottom=0.15)                                  
    #cax.cmap = self.colormap
    # if args.clog:
    #     cax.norm = matplotlib.colors.LogNorm(vmin=self.clim[0], vmax=self.clim[1])

    ax = plt.gca()
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_xlabel(xlabel, fontsize=14)
    cbar = fig.colorbar(cax)
    
    if args.part_selection_criterion != '' or args.tagfile_in != None:
        appstr += '_sel'

    raw.read_attrs()
    timestamp = '_%08.1f' % (raw.get_time()) 
    savename = savenamex \
                + '_' \
                + savenamey \
                + get_zeta_range_str(args.zeta_range) \
                + appstr \
                + timestamp

    saveas_png(fig, args.savepath, savename, verbose=True)
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser()
    r2h_subparsers = parser.add_subparsers(title="hist-type")

    r2h_parent_parser = raw2hist_parser()

    r2h_ssp = r2h_2d_subparser( subparsers=r2h_subparsers,
                                parent_parser=r2h_parent_parser)
    r2h_ssp.set_defaults(func=twoD)

    r2h_lsp = r2h_1d_subparser( subparsers=r2h_subparsers,
                                parent_parser=r2h_parent_parser)
    r2h_lsp.set_defaults(func=oneD)

    args = parser.parse_args()
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)


    if args.latexfont:
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif') 

    h5flist = H5FList(args.path, h5ftype='raw')
    flist = h5flist.get(stride=args.stride)

    for file in flist:
        raw = HiRAW(file)
        raw.read_data(verbose=args.verbose)

        if args.verbose:
            print( "The file contains %i particles." % raw.get_npart() ) 

        # Particle selection
        ts = TagSelect()
        ts.import_raw_obj(raw)
        if args.tagfile_in != None:
            # particle selection from tagfile_in
            ts.read_taglist(args.tagfile_in)
        else:
            # particle selection from criterion
            ts.select_part(args.part_selection_criterion)
        raw = ts.get_raw_select()
        print( "%i particles selected" % raw.get_npart() )
        if args.tagfile_out != None:
            ts.write_taglist(args.tagfile_out)

        if args.zeta_range != None:
            raw.select_zeta_range(args.zeta_range)
        if raw.get_npart() > 0:           
            args.func(raw, args)
        else:
            print('Error:\tNo particles in file (or range)!')
            sys.exit()                


if __name__ == "__main__":
    main()         