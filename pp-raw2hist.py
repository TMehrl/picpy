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
import pp_defs
from pp_h5dat import HiRAW
from pp_h5dat import H5FList
from pp_plt_tools import saveas_png
from pp_plt_tools import saveas_eps_pdf


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
    return parser


def r2h_1d_subparser(subparsers, parent_parser):
    parser = subparsers.add_parser( "1d", parents=[parent_parser],
                                    help="1D histogram")    
    # Line plot specific arguments
    parser.add_argument(  '--psv', '--phase-space-var',
                          action='store',
                          dest="psv",
                          metavar="VAR",
                          choices=[ 'x', 'y', 'z', 'px', 'py', 'pz'],
                          default='x',
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
                          help='Range of lineout.',
                          action='store',
                          dest="range",
                          metavar=('XMIN', 'XMAX'),
                          nargs=2,
                          type=float,
                          default=None)
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
                          metavar=('NBINS1', 'NBINS1'),
                          nargs=2,
                          type=int,
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


def oneD(raw, args):
    if args.psv == 'x':
       xlabel = r'$k_p x$' 
       psv = raw.x2
       savename = 'x'
    elif args.psv == 'z':
       xlabel = r'$k_p \zeta$' 
       psv = raw.x1

    if args.nbins == None:
        nbins = np.int(raw.get_npart()/1e2)
    else:
        nbins = args.nbins  

    hist, bin_edges = np.histogram(psv,bins=nbins,weights=raw.q)
    x_array = bin_edges[0:-1] + (bin_edges[1] - bin_edges[0])/2

    fig = plt.figure()
    plt.plot( x_array, hist )
    ax = plt.gca()    
    ax.set_xlabel(xlabel, fontsize=14)

    saveas_eps_pdf(fig, args.savepath, savename, h5plot=True, verbose=True, fformat='pdf')



def twoD(raw, args):
    # Do stuff  
    pass




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
    flist = h5flist.get()

    for file in flist:
        raw = HiRAW(file)
        raw.read_data()
        if args.zeta_range != None:
            raw.select_zeta_range(args.zeta_range)

        args.func(raw, args)


if __name__ == "__main__":
    main()         