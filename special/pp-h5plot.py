#!/usr/bin/env python3

import sys
import os
import math
import argparse
import numpy as np
import matplotlib
import itertools
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter

mypath = os.path.dirname(os.path.abspath( __file__ ))
incpath = os.path.split(mypath)[0] + '/inc'
sys.path.append(incpath)
from pp_h5dat import H5FList
from pp_h5dat import H5Plot
from pp_h5dat import mkdirs_if_nexist
from pp_plt_tools import saveas_eps_pdf,saveas_png
from pp_plt_tools import Colors

colors = Colors()

def h5plot_parser():

    desc = """This is the picpy postprocessing tool."""
    # Line vs-theo plot arguments
    parser = argparse.ArgumentParser( add_help=True,
                                      description=desc)
    parser.add_argument(  'paths',
                          metavar = 'PATHS',
                          nargs = '*',
                          help = 'Path(s) to hdf5 plot file.')
    parser.add_argument(  "-v", "--verbose",
                          dest = "verbose",
                          action="store_true",
                          default=True,
                          help = "Print info (Default).")
    parser.add_argument(  "-q", "--quiet",
                          dest = "verbose",
                          action = "store_false",
                          help = "Don't print info.")
    parser.add_argument(  "-s", "--save-path",
                          action="store",
                          dest="savepath",
                          metavar="PATH",
                          default=None,
                          help = """Path to which generated files will be saved.
                              (Default: %(default)s)""")
    parser.add_argument(  "-n", "--save-name",
                          action="store",
                          dest="savename",
                          metavar="NAME",
                          default=None,
                          help = """Name of saved file.
                              (Default: %(default)s)""")
    parser.add_argument(  "-f", "--format",
                          action='store',
                          dest="file_format",
                          metavar="FORMAT",
                          choices=[ 'pdf',
                                    'eps',
                                    'png'],
                          default='pdf',
                          help= """Format of output file (Default: %(default)s).""")
    parser.add_argument(  "--savenumpy",
                          dest = "savenumpy",
                          action = "store_true",
                          default=False,
                          help = "Saves x and y array of plots to numpy file (default: %(default)s).")
    parser.add_argument(  "--dpi",
                          action='store',
                          dest="dpi",
                          default=600,
                          type=int,
                          help= """Dots per inch for png output (default: %(default)s).""")
    parser.add_argument(  "-l", "--line-no",
                          action='store',
                          dest="line_no",
                          metavar="LINE-NUMBERS",
                          type=int,
                          default=None,
                          nargs='+',
                          help= """Selected line number for each provided file.""")
    parser.add_argument(  "--lab",
                          action='store',
                          dest="labels",
                          metavar="LABELS",
                          type=str,
                          default=None,
                          nargs='+',
                          help= """Labels for each selected line.""")
    parser.add_argument(  "--cno", "--color-number",
                          action='store',
                          dest="cno",
                          metavar="CNUMBERS",
                          type=int,
                          default=None,
                          nargs='+',
                          help= """Color number for each selected line.""")
    parser.add_argument(  "--cna", "--color-name",
                          action='store',
                          dest="cna",
                          metavar="CNAMES",
                          type=str,
                          default=None,
                          nargs='+',
                          help= """Color names for each line (choose amongst: %s)""" % colors.get_names_alpha_sorted())
    parser.add_argument(  "--show-palette", '--show-colors',
                          dest = "showpalette",
                          action="store_true",
                          default=False,
                          help = "Show color palette and exit.")
    parser.add_argument(  "--lstyle", "--line-style",
                          action='store',
                          dest="lstyle",
                          metavar="LINE-STYLE-NUMBER",
                          type=int,
                          default=None,
                          nargs='+',
                          help= """Line style number for each selected line. 0: '-', 1: '--', 2: '-.',  3: ':'""")
    parser.add_argument(  '--xlim',
                          help='x-range of plot.',
                          action='store',
                          dest="xlim",
                          metavar=('XMIN', 'XMAX'),
                          nargs=2,
                          type=float,
                          default=None)
    parser.add_argument(  '--ylim',
                          help='y-range of plot.',
                          action='store',
                          dest="ylim",
                          metavar=('YMIN', 'YMAX'),
                          nargs=2,
                          type=float,
                          default=None)
    parser.add_argument(  '--fig-size',
                          help='Size of figure in inch.',
                          action='store',
                          dest="figsize",
                          metavar=('WIDTH', 'HEIGHT'),
                          nargs=2,
                          type=float,
                          default=None)
    parser.add_argument(  '--xlab',
                          help='Set x label.',
                          action='store',
                          dest="xlab",
                          metavar='XLAB',
                          default=None)
    parser.add_argument(  '--ylab',
                          help='Set y label.',
                          action='store',
                          dest="ylab",
                          metavar='YLAB',
                          default=None)
    parser.add_argument(  "--xlog",
                          dest = "absxlog",
                          action="store_true",
                          default=False,
                          help = "Log x-axis (default: %(default)s).")
    parser.add_argument(  "--ylog",
                          dest = "absylog",
                          action="store_true",
                          default=False,
                          help = "Log y-axis (default: %(default)s).")
    parser.add_argument(  "--diff",
                          dest = "diff",
                          action="store_true",
                          default=False,
                          help = "Plot difference w.r.t. first line plot (default: %(default)s).")
    parser.add_argument(  "--abs",
                          dest = "abs",
                          action="store_true",
                          default=False,
                          help = "Plot absolute values of line plots (default: %(default)s).")
    parser.add_argument(  "--latexon",
                          dest = "latexon",
                          action="store_true",
                          default=False,
                          help = "Use LaTeX font (Default: %(default)s).")
    parser.add_argument(  "--xfac",
                          dest = "xfac",
                          action="store",
                          type=float,
                          default=1.0,
                          help = "Multiply y-values with specified factor (Default: %(default)s).")
    parser.add_argument(  "--yfac",
                          dest = "yfac",
                          action="store",
                          type=float,
                          default=1.0,
                          help = "Multiply y-values with specified factor (Default: %(default)s).")
    parser.add_argument(  "--noleg",
                          dest = "noleg",
                          action="store_true",
                          default=False,
                          help = "Disable legend (Default: %(default)s).")
    parser.add_argument(  "--flipleg",
                          dest = "flipleg",
                          action="store_true",
                          default=False,
                          help = "Flip legend entries (Default: %(default)s).")
    parser.add_argument(  "--legorder",
                          dest = "legorder",
                          type=int,
                          default=None,
                          nargs='+',
                          help= """Order of how the the legend entries should appear. E.g. '1 0 2' for three curves.""")
    parser.add_argument(  "--adj-left",
                          dest = "adjleft",
                          action="store",
                          type=float,
                          default=0.19,
                          help = "Adjust left plot margin (Default: %(default)s).")
    parser.add_argument(  "--adj-bottom",
                          dest = "adjbottom",
                          action="store",
                          type=float,
                          default=0.19,
                          help = "Adjust bottom plot margin (Default: %(default)s).")
    return parser



def main():

    parser = h5plot_parser()
    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    if args.showpalette:
        colors.show_palette()
        sys.exit(1)

    if args.figsize == None:
        figsize = (6,5)
    else:
        figsize = args.figsize

    fig = plt.figure(figsize=figsize)

    if (args.line_no != None) and (len(args.line_no) != len(args.paths)):
        print('ERROR: number of selected lines must be equal to the number of provided hdf5 files!')
        sys.exit(1)

    if (args.labels != None) and (len(args.labels) != len(args.paths)):
        print('ERROR: number of prepend labels must be equal to the number of provided hdf5 files!')
        sys.exit(1)

    h5lp = []

    linestyles = ['-', '--', '-.', ':']
    save_append_str = ''
    type_str = 'comp'

    if args.latexon:
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

    j = 0
    if args.diff:
        h5firstlp = H5Plot()
        h5firstlp.read(args.paths[0])
        type_str = 'diff'
        for (x, y, label, linestyle, color) in h5firstlp.get_line_plots():
            if j == 0:
                x0 = x
                y0 = y
            j += 1

    i = 0
    for path in args.paths:
        if not os.path.isfile(path):
            print('ERROR: File %s does not exist!' % path)
            sys.exit(1)

        h5lp.append(H5Plot())
        h5lp[i].read(path)

        if (args.cno == None) and (args.cna == None):
            argcolor = colors.get(i)
        elif (args.cno != None):
            argcolor = colors.get(args.cno[i])
        elif (args.cna != None):
            argcolor = colors.get(args.cna[i])

        j = 0
        for (x, y, label, linestyle, color) in h5lp[i].get_line_plots():

            # For normalizations etc.:
            x = args.xfac * x
            y = args.yfac * y

            if args.labels != None:
                save_append_str += '_' + args.labels[i]
                label = args.labels[i]

            if label == '_line0':
                label = path

            if args.latexon:
                if not '$' in label:
                    label = r'$\textrm{' + label + '}$'
                else:
                    label = r'%s' % label

            if args.diff:
                if (j == 0 and i == 0):
                    continue
                else:
                    y -= np.interp(x,x0,y0)

            if args.abs:
                y = np.abs(y)

            if args.lstyle == None:
                linestyle_code = linestyles[i%len(linestyles)]
            else:
                linestyle_code = linestyles[args.lstyle[i]]

            if (args.line_no != None):
                if args.line_no[i] == j:
                    plt.plot(x, y, label=label, linestyle=linestyle_code, color=argcolor)
            else:
                plt.plot(x, y, label=label, linestyle=linestyle_code, color=argcolor)
            j += 1
        i += 1

    ax = plt.gca()

    if args.xlim != None:
        ax.set_xlim(args.xlim[0],args.xlim[1])
    if args.ylim != None:
        ax.set_ylim(args.ylim[0],args.ylim[1])

    if args.absxlog:
        ax.set_xscale('log')
        ax.get_xaxis().get_major_formatter().labelOnlyBase = False

    if args.absylog:
        ax.set_yscale('log')
        ax.get_yaxis().get_major_formatter().labelOnlyBase = False

    if args.xlab != None:
        xlab = args.xlab
    else:
        xlab = h5lp[0].get_xlab()

    if args.ylab != None:
        ylab = args.ylab
    else:
        ylab = h5lp[0].get_ylab()

    if args.latexon:
        if xlab[0] != '$' or xlab[-1] != '$':
                    xlab = r'$' + xlab + '$'

    if args.latexon:
        if ylab[0] != '$' or ylab[-1] != '$':
                    ylab = r'$' + ylab + '$'

    ax.set_xlabel(xlab, fontsize=14)
    ax.set_ylabel(ylab, fontsize=14)
    handles, labels = ax.get_legend_handles_labels()
    if args.noleg != True:
        if args.flipleg == True:
            plt.legend(handles[::-1], labels[::-1],frameon=False)
            #plt.legend(flip(handles, 2), flip(labels, 2), ncol=2)
        elif args.legorder != None:
            handles_reordered = [handles[i] for i in args.legorder]
            labels_reordered = [labels[i] for i in args.legorder]
            plt.legend(handles_reordered, labels_reordered,frameon=False)
        else:
            plt.legend(frameon=False)
    plt.gcf().subplots_adjust(left=args.adjleft, bottom=args.adjbottom)
    if np.max(np.abs(y)) > 0.0:
        if not (-3.0 < math.log(np.max(np.abs(y)),10) < 3.0):
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
    plt.show()

    spath, fname  = os.path.split(args.paths[0])

    if args.savepath != None:
        spath = args.savepath

    if args.savename != None:
        save_name = args.savename
    else:
        save_name = fname + save_append_str + '_' + type_str

    if args.savenumpy == True:
        np.save(spath + save_name, [x,y])

    if args.file_format == 'png':
        saveas_png(fig, savepath=spath, savename=save_name, dpi=args.dpi)
    else:
        saveas_eps_pdf(fig, savepath=spath, savename=save_name, h5plot=True, verbose=True, fformat=args.file_format)
    plt.close(fig)

if __name__ == "__main__":
    main()
