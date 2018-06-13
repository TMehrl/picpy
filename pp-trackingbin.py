#!/usr/bin/env python3
import csv
import os
import glob
import numpy as np
import argparse
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import cm
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
    parser.add_argument(  "--2Dproj",
                          dest = "twodproj",
                          default = False,
                          help = "Plots a 2D projection of the particle" 
                          "trajectories instead of full 3D plot. \n"
                          "Use only if acquired data is based on a single slice! (default: %(default)s).")


    return parser


def plot_3D_colourline(x,y,z,c, minc, maxc):
    
    c = cm.jet((c-minc)/(maxc-minc))
    ax = plt.gca()
    for i in np.arange(len(x)-1):
        ax.plot([x[i],x[i+1]], [y[i],y[i+1]], [z[i],z[i+1]], c=c[i])
    
    return
    
def plot_2D_colourline(x,z,c, minc, maxc):
    
    c = cm.jet((c-minc)/(maxc-minc))
    ax = plt.gca()
    for i in np.arange(len(x)-1):
        ax.plot([x[i],x[i+1]], [z[i],z[i+1]], c=c[i])
    
    return


def main():

    parser = binSlab_parser() 
    args = parser.parse_args()

    # NHC=2
    proc_suffix_str = '_proc_'
    ppart_track_str = 'ppart_track'
    bin_fending = '.bin'
    
    for filepath in args.path:
        
        
        
        i_sl = np.int(filepath[(filepath.find(proc_suffix_str)+len(proc_suffix_str)):(len(filepath)-len(bin_fending))])
        #print(i_sl)

        if i_sl == 0:
            print('I AM ALIVE')
            filepath_base = filepath[0:(filepath.find(proc_suffix_str)+len(proc_suffix_str))]
            #print(filepath_base)
            Nsl = len(glob.glob(filepath_base + '*'))
            #print(Nsl)
        
            array = np.empty(0)
            for i in range(0,Nsl):
                filepath_slice = '%s%d%s' % (filepath_base, i, bin_fending )
                print('Reading: %s' % filepath_slice)
                array = np.append(array,np.fromfile(filepath_slice,dtype=np.float32))


            # array = np.empty(0)
            # 
            # filepath_slice = '%s%d%s' % (filepath_base, 1, bin_fending )
            # print('Reading: %s' % filepath_slice)
            # array = np.append(array,np.fromfile(filepath_slice,dtype=np.float32))
        
            # filepath_slice = '%s%d%s' % (filepath_base, 1, bin_fending )
            # print('Reading: %s' % filepath_slice)
            # array = np.append(array,np.fromfile(filepath_slice,dtype=np.float32))
            # filepath_base = filepath
            # array = np.empty(0)
            # for files in glob.glob(filepath_base + '/*'):
            #     #filepath_slice = '%s%d.bin' % (filepath_base, i )
            #     print('Reading: %s' % files)
            #     array = np.append(array,np.fromfile(files,dtype=np.float32))
        
            ''' reshape the binary data to a matrix with particle
            tag information in each row'''
            array = np.reshape(array, (int(len(array)/11),11))
            print(np.shape(array))
            ''' set up figure for 3D plotting '''
        
        
            ''' In order to simplify the splitting and plotting, the particle
            information matrix is sorted by the plasma species, the proc tag, then by the particle tag
            and finally by zeta '''
            indizes = np.lexsort((array[:,5], array[:,6], array[:,7], array[:,9])) 
            array = array[indizes]
        
            ''' split by different plasma species types '''
            w = np.split(array, np.where(np.diff(array[:, 9]) != 0)[0]+1) 
        
            for k in range(len(w)):
        
        
                fig = plt.figure()
                if args.twodproj:
                    ax = fig.add_subplot(111)
                else:
                    ax = fig.add_subplot(111, projection='3d')
                ''' get min and max value for universal colorbar later '''
                cmin = min(w[k][:,8])
                cmax = max(w[k][:,8])
        
                ''' Splitting the array into each particle trajectory by proc tag'''
                d= np.split(w[k], np.where(np.diff(w[k][:, 7]) != 0)[0]+1) 

                for i in range(len(d)):
        
                    ''' Splitting the array into each particle trajectory by part tag'''
                    e = np.split(d[i], np.where(np.diff(d[i][:, 6]) != 0)[0]+1) 
                    #print(np.shape(e))
                    #print('particle tags %i (always look at 1 and ignore 0) %i'% (i, e[16][1, 6]))
                    for j in range(len(e)):
                        x=e[j][:,0]
                        y=e[j][:,1]
                        z=e[j][:,5]
                        c=e[j][:,8]
                        # print(len(x))
                        # print(len(y))
                        print("laenge track proc tag %i part tag %i ist %i" %(i, j, len(z)))
                        
                        if args.twodproj:
                            plot_2D_colourline(z,x,c, cmin, cmax)
                        else:
                            plot_3D_colourline(z,y,x,c, cmin, cmax)
                        
                        #ax.plot(z, y, x, label='particle track')
        
                ''' Set colorbar ''' 
                norm = matplotlib.colors.Normalize(
                vmin=np.min(cmin),
                vmax=np.max(cmax))
                # choose a colormap
                c_m = cm.jet
                # create a ScalarMappable and initialize a data structure
                s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
                s_m.set_array([])
        
                cbar = plt.colorbar(s_m)
                cbar.ax.set_ylabel(r'$|u|$')
        
                ax.set_xlim(0, 300)
                ax.set_ylim(-6, 6)
                if not args.twodproj:
                    ax.set_zlabel(' x ')
                    ax.set_zlim(-6, 6)
                ax.set_xlabel(r'$\zeta$')
                ax.set_ylabel(' y ')
                
                #ax.legend()
        
                plt.show()
        #     # fig = plt.figure()
        #     # cax = plt.plot( array)
        #     # 
        #     # ax = plt.gca()
        #     # ax.set_ylabel(r'$\eta (\zeta)$', fontsize=14)
        #     # ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
        #     # 
        #     # if args.ifshow:
        #     #   plt.show()
        #     # 
        #     # head, savename = os.path.split(filepath)
        #     # 
        #     # mkdirs_if_nexist(args.savepath)
        #     # 
        #     # print('Saving: ' + args.savepath + '/' + savename + '.eps')
        #     # 
        #     # fig.savefig( args.savepath + '/' + savename + '.eps',
        #     #         format='eps')    

if __name__ == "__main__":
    main()  