#!/usr/bin/env python3
import csv
import os
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib import cm
import pp_defs
from pp_h5dat import mkdirs_if_nexist

def binSlab_parser():

    desc = """This is the picpy postprocessing tool."""

    parser = argparse.ArgumentParser( description=desc)
    parser.add_argument(  'path',
                          metavar = 'PATH',
                          nargs = '*',
                          help = 'Path to slice files.')
    parser.add_argument(  '-o', '--outp-sum-file',
                          action='store',
                          metavar = 'OUTPUT_SUM',
                          dest = "outsumfile",
                          required=True,
                          help = 'Path to output summary file.')    
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
    parser.add_argument(  "-c", "--code",
                          action="store",
                          dest="piccode",
                          choices = [pp_defs.code.hipace, pp_defs.code.osiris,],
                          default = pp_defs.code.hipace,
                          help="PIC code (default: %(default)s).")
    parser.add_argument(  "--save-path",
                          action="store",
                          dest="savepath",
                          metavar="SPATH",
                          default='./plots_debug/',
                          help = """Path to which generated files will be saved.
                              (default: %(default)s)""")
    parser.add_argument(  "-p", "--partition",
                          help='Virtual topology partitioning',
                          action='store',
                          dest="partition",
                          metavar=('NPROCX0', 'NPROCX1', 'NPROCX2'),
                          nargs=3,
                          type=int,
                          default=(1,1,1))
    parser.add_argument(  "--nohalo",
                          dest = "nohalo",
                          action="store_true",
                          default=False,
                          help = "If data is w/o halo (default: %(default)s).")                                                       
    return parser


def main():

    parser = binSlab_parser() 
    args = parser.parse_args()

    NHC=2

    #cblims = [-1e-1 1e-1];

    grid_info = np.zeros((3, 3))
    
    with open(args.outsumfile, 'r') as csvfile:
        csvread = csv.reader(csvfile, delimiter='\t')
        i=0
        for line in csvread:
            if (i<3):
                for j in range(len(line)):
                    grid_info[i,j] = float(line[j])
            i += 1        

    print(grid_info)

    Nx_wo_halo = int(grid_info[1,2])
    Ny_wo_halo = int(grid_info[2,2])

    dx=(grid_info[1,1] - grid_info[1,0]) / (Nx_wo_halo - 1)
    dy=(grid_info[2,1] - grid_info[2,0]) / (Ny_wo_halo - 1)

    Nx_loc=np.int(Nx_wo_halo/args.partition[1])
    Ny_loc=np.int(Ny_wo_halo/args.partition[2])

    Nx_loc_incl_halo = Nx_loc + 2*NHC;
    Ny_loc_incl_halo = Ny_loc + 2*NHC;

    x_array_wo_halo=np.linspace(grid_info[1,0],grid_info[1,1],Nx_loc)
    y_array_wo_halo=np.linspace(grid_info[2,0],grid_info[2,1],Ny_loc)
    
    x_array_incl_halo = np.linspace(grid_info[1,0] - NHC*dx, 
                                    grid_info[1,1] + NHC*dx,
                                    Nx_loc_incl_halo)
    y_array_incl_halo = np.linspace(grid_info[2,0] - NHC*dy,
                                    grid_info[2,1] + NHC*dy,
                                    Ny_loc_incl_halo)
    

    for filepath in args.path:

        if args.nohalo:
            Nx=Nx_loc
            Ny=Ny_loc

            x_array=x_array_wo_halo
            y_array=y_array_wo_halo

        else:
            Nx=Nx_loc_incl_halo
            Ny=Ny_loc_incl_halo

            x_array=x_array_incl_halo
            y_array=y_array_incl_halo

        M_1D = np.fromfile(filepath,dtype=np.float32)
        M = np.transpose(M_1D.reshape((Nx, Ny)))

        fig = plt.figure()
        cax = plt.pcolormesh( x_array,
                              y_array,
                              M,
                              cmap=cm.PuBu)

        ax = plt.gca()
        ax.set_ylabel(r'$k_p y$', fontsize=14)
        ax.set_xlabel(r'$k_p x$', fontsize=14)
        cbar = fig.colorbar(cax)
        #cbar.ax.set_ylabel( gen_pretty_grid_name( self.g3d.name ), fontsize=14 )
        if args.ifshow:
            plt.show()

        head, savename = os.path.split(filepath)

        mkdirs_if_nexist(args.savepath)    
        fig.savefig( args.savepath + '/' + savename + '.png',
                  format='png',
                  dpi=600)    

if __name__ == "__main__":
    main()  