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
from pp_plt_tools import saveas_png


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
                          default=None,
                          help = """Path to which generated files will be saved.
                              (Default: %(default)s)""")
    parser.add_argument(  "-f", "--format",
                          action='store',
                          dest="file_format",
                          metavar="FORMAT",
                          choices=[ 'png',
                                    'pdf',
                                    'eps',],
                          default='png',
                          help= """Format of output file (Default: %(default)s).""")
    parser.add_argument(  "-l", "--line-no",
                          action='store',
                          dest="line_no",
                          metavar="LINE-NUMBERS",
                          type=int,                      
                          default=None,
                          nargs='+',
                          help= """Selected line number for each provided file.""")
    parser.add_argument(  "--data2",
                        action='store',
                        dest="data2",
                        metavar="DATA2",
                        nargs = '*',
                        default=None,
                        help= """Path to the second data set""")
    parser.add_argument(  "--lab",
                          action='store',
                          dest="labels",
                          metavar=('LABEL1', 'LABEL2'),
                          type=str,                      
                          default=None,
                          nargs=2,
                          help= """Labels for the two lines.""")
    parser.add_argument(  "--cno", "--color-number",
                          action='store',
                          dest="cno",
                          metavar=("CNUMBER1","CNUMBER2"),
                          type=int,                      
                          default=None,
                          nargs=2,
                          help= """Color number for each selected line.""")
    parser.add_argument(  "--lstyle", "--line-style",
                          action='store',
                          dest="lstyle",
                          metavar=("LINE-STYLE-NUMBER1", "LINE-STYLE-NUMBER2" ),
                          type=int,                      
                          default=None,
                          nargs=2,
                          help= """Line style number for both lines. 0: '-', 1: '--', 2: '-.',  3: ':'""")
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
    parser.add_argument(  '--y2lim',
                          help='y-range of plot 2.',
                          action='store',
                          dest="y2lim",
                          metavar=('YMIN', 'YMAX'),
                          nargs=2,
                          type=float,
                          default=None)  
    parser.add_argument(  '--xlab',
                          help='Set x label.',
                          action='store',
                          dest="xlab",
                          metavar='XLAB1',
                          default=None)  
    parser.add_argument(  '--ylab',
                          help='Set y label.',
                          action='store',
                          dest="ylab",
                          nargs=2,
                          metavar=('YLAB1','YLAB2'),
                          default=None)                                                         
    parser.add_argument(  "--ylog",
                          dest = "absylog",
                          action="store_true",
                          default=False,
                          help = "Plot abs log of y-data (default: %(default)s).")
    parser.add_argument(  "--y2log",
                        dest = "absy2log",
                        action="store_true",
                        default=False,
                        help = "Plot abs log of y-data2 (default: %(default)s).") 
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
    parser.add_argument(  "--dpi",
                          action='store',
                          dest="dpi",
                          default=300,
                          type=int,
                          help= """Dots per inch for png output (default: %(default)s).""")
    parser.add_argument(  "--n0",
                          action='store',
                          dest="n0",
                          default=1e23,
                          type=float,
                          help= """Plasma density n0 in order to calculate length of simulation (default: %(default)s).""")
    parser.add_argument(  "--maketitle",
                          dest = "maketitle",
                          action="store_true",
                          default=False,
                          help = "Writes title with time and length (Default: %(default)s).") 
    return parser

 

def main():


    parser = h5plot_parser()
    args = parser.parse_args()

    if args.maketitle:
        ''' Calculating plasma frequency for length '''
        ELECTRON_CHARGE_IN_COUL   =    1.60217657e-19
        ELECTRON_MASS_IN_KG       =    9.10938291e-31
        VAC_PERMIT_FARAD_PER_M    =    8.854187817e-12
        SPEED_OF_LIGHT_IN_M_PER_S =    299792458.0

        omega_p = np.sqrt( args.n0 * (ELECTRON_CHARGE_IN_COUL**2)/ (VAC_PERMIT_FARAD_PER_M * ELECTRON_MASS_IN_KG));
        skin_depth = SPEED_OF_LIGHT_IN_M_PER_S/omega_p
        
        

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    saveformat = args.file_format

    h5lp = []

    linestyles = ['-', '--', '-.', ':']
    save_append_str = ''
    type_str = 'comp'

    if args.latexon:
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif') 

    j = 0

    i = 0
    
    
    
    for path in args.paths:
        fig = plt.figure(figsize=(7,5))
        if not os.path.isfile(path):
            print('ERROR: File %s does not exist!' % path)
            sys.exit(1)          
            
        h5lp.append(H5Plot())
        h5lp[i].read(path)
        timestamp = path.split("_")[-1].split(".h5")[0]
        
        if args.data2:
            for files in args.data2:
                if timestamp in files:
                    h5secondlp = H5Plot()
                    h5secondlp.read(files)

        
        
        
        if args.cno == None:
            argcolor = plt.cm.tab20(6)
        else:
            argcolor = plt.cm.tab20(args.cno[0])
        
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


            if args.abs:
                y = np.abs(y)

            if args.lstyle == None:
                linestyle_code = linestyles[0]
            else:
                linestyle_code = linestyles[args.lstyle[0]]

            
            if args.absylog:
                plt.semilogy( x, y, label=label, linestyle=linestyle_code, color=argcolor, zorder=20)
            else:    
                #plt.plot(x, y, label=label, linestyle=linestyle, color=color)
                plt.plot(x, y, label=label, linestyle=linestyle_code, color=argcolor, zorder=20)      
            

        ax = plt.gca()
        ax.tick_params('y', colors=argcolor)
        plt.gcf().subplots_adjust(left=0.15, bottom=0.15, right=1-0.15)       
        if not (-3.0 < math.log(np.max(abs(y)),10) < 3.0):
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
            plt.gcf().subplots_adjust(left=0.18) 
        ax2 = ax.twinx()
        
        ''' PLOTTING THE SECOND PLOT '''
        if args.cno == None:
            argcolor = plt.cm.tab20(1)
        else:
            argcolor = plt.cm.tab20(args.cno[1])
        ax2.tick_params('y', colors=argcolor)
        j = 0
        for (x2, y2, label, linestyle, color) in h5secondlp.get_line_plots():
            if args.labels != None:
                save_append_str += '_' + args.labels[i]
                label = args.labels[i]

            if label == '_line0':
                label = path

            if args.latexon:
                if label[0] != '$' or label[-1] != '$':
                    label = r'$\textrm{' + label + '}$'

            if args.abs:
                y2 = np.abs(y2)

            if args.lstyle == None:
                linestyle_code = linestyles[0]
            else:
                linestyle_code = linestyles[args.lstyle[1]]

            
            if args.absy2log:
                plt.semilogy( x2, y2, label=label, linestyle=linestyle_code, color=argcolor, zorder=0)
            else:    
                #plt.plot(x, y, label=label, linestyle=linestyle, color=color)
                plt.plot(x2, y2, label=label, linestyle=linestyle_code, color=argcolor, zorder=0)      
            



        if args.xlim != None:
            ax.set_xlim(args.xlim[0],args.xlim[1])
        if args.ylim != None:
            ax.set_ylim(args.ylim[0],args.ylim[1])
        if args.y2lim != None:
            ax2.set_ylim(args.y2lim[0],args.y2lim[1])
            
        if args.xlab != None:
            xlab = args.xlab
        else:
            xlab = h5lp[0].get_xlab()
    

        if args.ylab != None:
            ylab = args.ylab
        else:
            ylab = []
            ylab.append(h5lp[0].get_ylab())
            ylab.append(h5secondlp.get_ylab())
            
        if args.latexon:
            if xlab[0] != '$' or xlab[-1] != '$':
                        xlab = r'$' + xlab + '$'
                        
                        
        if args.latexon:
            if ylab[0,0] != '$' or ylab[0,-1] != '$':
                        ylab[0] = r'$' + ylab[0] + '$'
                        ylab[1] = r'$' + ylab[1] + '$'

        ax.set_xlabel(xlab, fontsize=14)
        ax.set_ylabel(ylab[0], fontsize=14)
        ax2.set_ylabel(ylab[1], fontsize=14)   
        
        
        ax.set_zorder(ax2.get_zorder()+1) # put ax in front of ax2
        ax.patch.set_visible(False) # hide the 'canvas' 
        # handles, labels = ax.get_legend_handles_labels()
        # #plt.legend(flip(handles, 2), flip(labels, 2), ncol=2)
        # plt.legend(frameon=False)
        #plt.gcf().subplots_adjust(left=0.15, bottom=0.15, right=0.15)       
        if not (-3.0 < math.log(np.max(abs(y2)),10) < 3.0):
            ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
            plt.gcf().subplots_adjust(right=1-0.18)          
        #plt.show()

        ''' make the title '''
        if args.maketitle:
            time = np.around(float(timestamp))
            distance = time*skin_depth*1e2
            title = 'Time: ' + str(time) + r'$\,\omega_p^{-1}$ ' + '   Distance: ' + str(np.around(distance,2)) + r'$\,$cm'
            plt.title(title)
            

        spath, fname  = os.path.split(args.paths[0])
        if args.savepath != None:
            spath = args.savepath

        save_name = type_str + '_' + fname + save_append_str + '_' +  timestamp
        
        
        if saveformat==args.file_format:
            saveas_png(fig, spath, save_name, args.dpi)
        else:
            saveas_eps_pdf(fig, savepath=spath, savename=save_name)
        
        plt.close(fig)

if __name__ == "__main__":
    main()
