#!/usr/bin/env python3

import sys
import os
import math
import argparse
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import itertools
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter
from pp_h5dat import H5FList
from pp_h5dat import H5Plot
from pp_h5dat import mkdirs_if_nexist



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
    linestyles = ['-', '--', '-.', ':']
    save_append_str = ''
    
  
    for path in args.paths:
        if not os.path.isfile(path):
            print('ERROR: File %s does not exist!' % path)
            sys.exit(1)             
        
        h5lp.append(H5Plot())
        h5lp[i].read(path)
        j = 0 
        for (x, y, label, linestyle, color) in h5lp[i].get_line_plots():
            if args.labels != None:
                save_append_str += '_' + args.labels[i]
                label = args.labels[i]
            
            if label == '_line0':
                label = path

            if (args.line_no != None):
                if args.line_no[i] == j:
                    if args.absylog:
                        plt.semilogy( x, y, label=label, linestyle=linestyles[i%len(linestyles)])                        
                    else:
                        plt.plot(x, y, label=label, linestyle=linestyles[i%len(linestyles)])
            else:
                if args.absylog:
                    plt.semilogy( x, y, label=label, linestyle=linestyles[i%len(linestyles)])
                else:    
                    #plt.plot(x, y, label=label, linestyle=linestyle, color=color)
                    plt.plot(x, y, label=label, linestyle=linestyles[i%len(linestyles)])      
            j += 1
        i += 1
    
    ax = plt.gca()
    ax.set_xlabel(h5lp[0].get_xlab(), fontsize=14)
    ax.set_ylabel(h5lp[0].get_ylab(), fontsize=14)    
    handles, labels = ax.get_legend_handles_labels()
    #plt.legend(flip(handles, 2), flip(labels, 2), ncol=2)
    #plt.legend()

    # if not (-3.0 < math.log(np.max(abs(y)),10) < 3.0):
    #     ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
    #     plt.gcf().subplots_adjust(left=0.18)
    # else:
    #     plt.gcf().subplots_adjust(left=0.15)              
    plt.show()
    fname, fext = os.path.splitext(args.paths[0])
    save_path_name = fname + save_append_str + '_comp' + '.' + args.file_format
    fig.savefig(save_path_name, format=args.file_format)
    plt.close(fig)
    if args.verbose: 
        sys.stdout.write('Saved: %s\n' % save_path_name)
        sys.stdout.flush()

    if args.diff:
        print('DIFFERENCE PLOTS: TO BE IMPLEMENTED!!!')


if __name__ == "__main__":
    main()
