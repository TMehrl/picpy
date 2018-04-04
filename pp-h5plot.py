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
                          default=parsedefs.save.path + '/g3d-line-vs-theo',
                          help = """Path to which generated files will be saved.
                              (Default: %(default)s)""")
    parser.add_argument(  "-f", "--format",
                          action='store',
                          dest="file_format",
                          metavar="FORMAT",
                          choices=[ 'png',
                                    'pdf',
                                    'eps',],
                          default='eps',
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
                          help= """Labels to be prepended for each provided file.""")    
    return parser

 

def main():

    parser = h5plot_parser()
    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    fig = plt.figure()

    if (args.line_no != None) and (len(args.line_no) != len(args.paths)):
        print('ERROR: number of selected lines must be equal to the number of provided hdf5 files!')
        sys.exit(1)        

    if (args.labels != None) and (len(args.labels) != len(args.paths)):
        print('ERROR: number of prepend labels must be equal to the number of provided hdf5 files!')
        sys.exit(1)  
    
    h5lp = []
    i = 0
    for path in args.paths:
        h5lp.append(H5Plot())
        h5lp[i].read(path)
        j = 0   
        for (x, y, label, linestyle, color) in h5lp[i].get_line_plots():
            if args.labels != None:
                label = args.labels[i] + ' ' + label          
            if (args.line_no != None):
                if args.line_no[i] == j:
                    plt.plot(x, y, label=label, linestyle=linestyle)
            else:
                plt.plot(x, y, label=label, linestyle=linestyle, color=color)      
            j += 1
        i += 1

    ax = plt.gca()
    ax.set_ylabel(h5lp[0].get_xlab(), fontsize=14)
    ax.set_xlabel(h5lp[0].get_ylab(), fontsize=14)    
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(flip(handles, 2), flip(labels, 2), ncol=2)           
    plt.show()
    fname, fext = os.path.splitext(args.paths[0])
    fig.savefig(  fname + '_h5plot.' + args.file_format ,
                      format=args.file_format)
    plt.close(fig)

if __name__ == "__main__":
    main()
