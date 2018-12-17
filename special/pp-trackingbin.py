#!/usr/bin/env python3
import csv
import os
import glob
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib import cm

mypath = os.path.dirname(os.path.abspath( __file__ ))
incpath = os.path.split(mypath)[0] + '/inc'
sys.path.append(incpath)
import pp_defs
from pp_h5dat import mkdirs_if_nexist
from mpl_toolkits.mplot3d import Axes3D


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
                          required=False,
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
    parser.add_argument(  "--nohalo",
                          dest = "nohalo",
                          action="store_true",
                          default=True,
                          help = "If data is w/o halo (default: %(default)s).")                                                       
    return parser


def main():

    parser = binSlab_parser() 
    args = parser.parse_args()

    # NHC=2
    slice_suffix_str = '_izeta_'


    # grid_info = np.zeros((3, 3))
    
    # with open(args.outsumfile, 'r') as csvfile:
    #     csvread = csv.reader(csvfile, delimiter='\t')
    #     i=0
    #     for line in csvread:
    #         if (i<3):
    #             for j in range(len(line)):
    #                 grid_info[i,j] = float(line[j])
    #         i += 1        

    # print(grid_info)

    # Nx = int(grid_info[1,2])
    # Ny = int(grid_info[2,2])

    # dx=(grid_info[1,1] - grid_info[1,0]) / (Nx - 1) # Check this!!!!
    # dy=(grid_info[2,1] - grid_info[2,0]) / (Ny - 1) # Check this!!!!

    # Nx_incl_halo = Nx + 2*NHC;
    # Ny_incl_halo = Ny + 2*NHC;


    for filepath in args.path:

        # if args.nohalo:
        #     nhc = 0
        # else:
        #     nhc = NHC

        i_sl = 0 #np.int(filepath[(filepath.find(slice_suffix_str)+len(slice_suffix_str)):len(filepath)])
        print(i_sl)
        # IDx = np.int(np.floor(slID/args.partition[2]))
        # IDy = slID - args.partition[2] * IDx

        # ix_start = np.int((Nx * IDx/args.partition[1]))
        # ix_end = np.int((Nx * (IDx+1)/args.partition[1]))

        # iy_start = np.int((Ny * IDy/args.partition[2]))
        # iy_end = np.int((Ny * (IDy+1)/args.partition[2]))

        # my_Nx=np.int(Nx/args.partition[1] + 2*nhc )
        # my_Ny=np.int(Ny/args.partition[2] + 2*nhc )
    
        # my_x_array = np.linspace( grid_info[1,0] + ix_start*dx - nhc*dx, 
        #                           grid_info[1,0] + ix_end*dx + nhc*dx,
        #                           my_Nx)
        # my_y_array = np.linspace( grid_info[2,0] + iy_start*dy - nhc*dy,
        #                           grid_info[2,0] + iy_end*dy + nhc*dy,
        #                           my_Ny)

        if i_sl == 0:

            filepath_base = filepath[0:(filepath.find(slice_suffix_str)+len(slice_suffix_str))]
            Nsl = len(glob.glob(filepath_base + '*'))

            array = np.empty(0)
            for i in range(0,Nsl):
                filepath_slice = '%s%04d.bin' % (filepath_base, i )
                print('Reading: %s' % filepath_slice)
                array = np.append(array,np.fromfile(filepath_slice,dtype=np.float32))
            
            array = np.reshape(array, (int(len(array)/7),7))
            print(array)
            print(np.shape(array))
            
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            x=array[:,0]
            y=array[:,1]
            z=np.arange(len(array))
            print(len(x))
            print(len(y))
            print(len(z))
            ax.plot(x, z, label='x')
            ax.plot(y,z, label='y')
            ax.legend()
            
            plt.show()
            # fig = plt.figure()
            # cax = plt.plot( array)
            # 
            # ax = plt.gca()
            # ax.set_ylabel(r'$\eta (\zeta)$', fontsize=14)
            # ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
            # 
            # if args.ifshow:
            #   plt.show()
            # 
            # head, savename = os.path.split(filepath)
            # 
            # mkdirs_if_nexist(args.savepath)
            # 
            # print('Saving: ' + args.savepath + '/' + savename + '.eps')
            # 
            # fig.savefig( args.savepath + '/' + savename + '.eps',
            #         format='eps')    

if __name__ == "__main__":
    main()  