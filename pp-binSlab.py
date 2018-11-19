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
    slID_str = 'slID_'

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

    Nx = int(grid_info[1,2])
    Ny = int(grid_info[2,2])

    dx=(grid_info[1,1] - grid_info[1,0]) / (Nx - 1) # Check this!!!!
    dy=(grid_info[2,1] - grid_info[2,0]) / (Ny - 1) # Check this!!!!

    Nx_incl_halo = Nx + 2*NHC;
    Ny_incl_halo = Ny + 2*NHC;


    
    for filepath in args.path:

        if args.nohalo:
            nhc = 0
        else:
            nhc = NHC

        slID = np.int(filepath[(filepath.find(slID_str)+len(slID_str)):len(filepath)])

        IDx = np.int(np.floor(slID/args.partition[2]))
        IDy = slID - args.partition[2] * IDx

        ix_start = np.int((Nx * IDx/args.partition[1]))
        ix_end = np.int((Nx * (IDx+1)/args.partition[1]))

        iy_start = np.int((Ny * IDy/args.partition[2]))
        iy_end = np.int((Ny * (IDy+1)/args.partition[2]))

        my_Nx=np.int(Nx/args.partition[1] + 2*nhc )
        my_Ny=np.int(Ny/args.partition[2] + 2*nhc )
    
        my_x_array = np.linspace( grid_info[1,0] + ix_start*dx - nhc*dx, 
                                  grid_info[1,0] + ix_end*dx + nhc*dx,
                                  my_Nx)
        my_y_array = np.linspace( grid_info[2,0] + iy_start*dy - nhc*dy,
                                  grid_info[2,0] + iy_end*dy + nhc*dy,
                                  my_Ny)

        array_flat = np.fromfile(filepath,dtype=np.float32)
        M = np.transpose(array_flat.reshape((my_Nx, my_Ny)))

        cblim = [0.0, 0.0]

        if np.amin(M) < 0 and np.amax(M) > 0:
            cblim[0] = -max(abs(np.amin(M)),np.amax(M))
            cblim[1] = max(abs(np.amin(M)),np.amax(M))
            colormap = cm.BrBG
        elif np.amin(M) < 0 and np.amax(M) <= 0:
            cblim[0] = np.amin(M)
            cblim[1] = 0
            colormap = cm.gist_yarg
        elif np.amin(M) >= 0 and np.amax(M) > 0:
            cblim[0] = 0
            cblim[1] = np.amax(M)
            colormap = cm.gist_gray
        else:
            cblim[0] = np.amin(M)
            cblim[1] = np.amax(M)
            colormap = cm.PuBu          

        fig = plt.figure()
        cax = plt.pcolormesh( my_x_array,
                              my_y_array,
                              M,
                              vmin=cblim[0], vmax=cblim[1],                              
                              cmap=colormap)

        ax = plt.gca()
        ax.set_ylabel(r'$k_p y$', fontsize=14)
        ax.set_xlabel(r'$k_p x$', fontsize=14)
        ax.set_aspect('equal')
        cbar = fig.colorbar(cax)

        if args.ifshow:
            plt.show()

        head, savename = os.path.split(filepath)

        mkdirs_if_nexist(args.savepath)    
        fig.savefig( args.savepath + '/' + savename + '.png',
                  format='png',
                  dpi=600)    

if __name__ == "__main__":
    main()  