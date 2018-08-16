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
from pp_h5dat import H5FList
from pp_h5dat import H5Plot
from pp_h5dat import mkdirs_if_nexist
from pp_plt_tools import saveas_eps_pdf


def flip(items, ncol):
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])


# class PositionalAction(argparse.Action):
#     def __call__(self,parser,namespace,values,option_string=None):
#         lst = getattr(namespace,self.dest)
#         lst.append(values)
#         parser.last_positional_values = lst
#         all_positional = getattr(namespace,'all_positional',[])
#         all_positional.append(lst)
#         namespace.all_positional = all_positional

# class AssociateAction(argparse.Action):
#     def __call__(self,parser,namespace,values,option_string=None):
#         try:
#             parser.last_positional_values.append(values)
#         except AttributeError:
#             pass

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
                          default='plots/g3d-line-vs-theo',
                          help = """Path to which generated files will be saved.
                              (Default: %(default)s)""")
    parser.add_argument(  "-f", "--format",
                          action='store',
                          dest="file_format",
                          metavar="FORMAT",
                          choices=[ 'png',
                                    'pdf',
                                    'eps',],
                          default='pdf',
                          help= """Format of output file (Default: %(default)s).""")
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
    parser.add_argument(  "--lstyle", "--line-style",
                          action='store',
                          dest="lstyle",
                          metavar="LINE-STYLE-NUMBER",
                          type=int,                      
                          default=None,
                          nargs='+',
                          help= """Line style number for each selected line. 0: '-', 1: '--', 2: '-.',  3: ':'""")
    # parser.add_argument(  "--lstyle", "--line-style",
    #                       action='store',
    #                       dest="lstyle",
    #                       metavar="LINE-STYLE",                          
    #                       choices=[ '-', '-.', ':'],
    #                       nargs='+',
    #                       default=None,
    #                       help= """Line style number for each selected line (default: %(default)s).""")
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
    parser.add_argument(  "--ylog",
                          dest = "absylog",
                          action="store_true",
                          default=False,
                          help = "Plot abs log of y-data (default: %(default)s).") 
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
    return parser

 

def main():

    parser = h5plot_parser()
    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    fig = plt.figure(figsize=(6,5))

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

        if args.cno == None:
            argcolor = plt.cm.tab20c(i)
        else:
            argcolor = plt.cm.tab20c(args.cno[i])
        
        j = 0
        for (x, y, label, linestyle, color) in h5lp[i].get_line_plots():
            if args.labels != None:
                save_append_str += '_' + args.labels[i]
                label = args.labels[i]

            if label == '_line0':
                label = path

            if args.latexon:
                if label[0] != '$' or label[-1] != '$':
                    label = r'$\textrm{' + label + '}$'

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
                    if args.absylog:
                        plt.semilogy( x, y, label=label, linestyle=linestyle_code, color=argcolor)                        
                    else:
                        plt.plot(x, y, label=label, linestyle=linestyle_code, color=argcolor)
            else:
                if args.absylog:
                    plt.semilogy( x, y, label=label, linestyle=linestyle_code, color=argcolor)
                else:    
                    #plt.plot(x, y, label=label, linestyle=linestyle, color=color)
                    plt.plot(x, y, label=label, linestyle=linestyle_code, color=argcolor)      
            j += 1
        i += 1

    ax = plt.gca()

    if args.xlim != None:
        ax.set_xlim(args.xlim[0],args.xlim[1])
    if args.ylim != None:
        ax.set_ylim(args.ylim[0],args.ylim[1])

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
    #plt.legend(flip(handles, 2), flip(labels, 2), ncol=2)
    plt.legend(frameon=False)
    plt.gcf().subplots_adjust(left=0.15, bottom=0.15)       
    if not (-3.0 < math.log(np.max(abs(y)),10) < 3.0):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)          
    plt.show()
    spath, fname  = os.path.split(args.paths[0])

    save_name = fname + save_append_str + '_' + type_str
    saveas_eps_pdf(fig, savepath=spath, savename=save_name)
    plt.close(fig)

if __name__ == "__main__":
    main()
