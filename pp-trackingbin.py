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

from matplotlib.colors import LinearSegmentedColormap

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

    parser.add_argument(  "-t", "--track_color",
                          dest = "track_color",
                          default = "u_tot",
                          help = "Select the attribute which is displayed as the color of the track \n"
                                 "Default is the total momentum u_tot, other available options are: \n"
                                 "beta_z")



    return parser


def display_cmap(cmap):
    plt.imshow(np.linspace(0, 100, 256)[None, :],  aspect=25,    interpolation='nearest', cmap=cmap) 
    plt.axis('off')
    plt.show()


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

def plot_3D_colourline_beta(x,y,z,c):
    basic_cols=['#2020ff' ,'#808080', '#ff2020']
    basic_cols=['#0000ff' , '#00ffff','#808080', '#ffff00', '#ff0000']
    my_cmap = LinearSegmentedColormap.from_list('mycmap', basic_cols)
    c = my_cmap((c+1)/(2))
    ax = plt.gca()
    for i in np.arange(len(x)-1):
        ax.plot([x[i],x[i+1]], [y[i],y[i+1]], [z[i],z[i+1]], c=c[i])
    
    return
    
def plot_2D_colourline_beta(x,z,c):
    basic_cols=['#2020ff' ,'#808080', '#ff2020']
    basic_cols=['#0000ff' , '#00ffff','#808080', '#ffff00', '#ff0000']
    my_cmap = LinearSegmentedColormap.from_list('mycmap', basic_cols)
    c = my_cmap((c+1)/(2))
    ax = plt.gca()
    for i in np.arange(len(x)-1):
        ax.plot([x[i],x[i+1]], [z[i],z[i+1]], c=c[i])
    
    return

def main():


    basic_cols=['#75b765', '#808080', '#ffd700']
    basic_cols=['#2020ff' ,'#808080', '#ff2020']
    basic_cols=['#0000ff' , '#00ffff','#808080', '#ffff00', '#ff0000']
    my_cmap=LinearSegmentedColormap.from_list('mycmap', basic_cols)
    # display_cmap(my_cmap)
    
    parser = binSlab_parser() 
    args = parser.parse_args()

    # NHC=2
    proc_suffix_str = '_proc_'
    ppart_track_str = 'ppart_track'
    bin_fending = '.bin'
    
    for filepath in args.path:

        filename = ''
        slash = '/'
        if slash in filepath:
            dirname, filename = filepath.split('/')
        else:
            dirname = filepath
        
        array = np.empty(0)
        for files in os.listdir(dirname):
            if filename != '':
                if filename in files:
                    print('Reading: %s/%s' % (dirname, files))
                    filepath = dirname + slash + files
                    array = np.append(array,np.fromfile(filepath,dtype=np.float32))
            else:
                print('Reading: %s/%s' % (dirname, files))
                filepath = dirname + slash + files
                array = np.append(array,np.fromfile(filepath,dtype=np.float32))
        
        # i_sl = np.int(filepath[(filepath.find(proc_suffix_str)+len(proc_suffix_str)):(len(filepath)-len(bin_fending))])
        # #print(i_sl)
        # i_sl =0
        # if i_sl == 0:
        #     print('I AM ALIVE')
        #     filepath_base = filepath[0:(filepath.find(proc_suffix_str)+len(proc_suffix_str))]
        #     #print(filepath_base)
        #     Nsl = len(glob.glob(filepath_base + '*'))
        #     #print(Nsl)
        # 
        #     array = np.empty(0)
        #     for i in range(0,Nsl):
        #         filepath_slice = '%s%d%s' % (filepath_base, i, bin_fending )
        #         print('Reading: %s' % filepath_slice)
        #         array = np.append(array,np.fromfile(filepath_slice,dtype=np.float32))
        # 

            # array = np.empty(0)
            # 
            # filepath_slice = '%s%d%s' % (filepath_base, 2, bin_fending )
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
            if args.track_color == "u_tot":
                cmin = min(w[k][:,8])
                cmax = max(w[k][:,8])
                chosencmap = cm.jet
            elif args.track_color == "beta_z":
                cmin = -1
                cmax = 1
                chosencmap = my_cmap
            elif args.track_color == "beta_y":
                cmin = -1
                cmax = 1
                chosencmap = my_cmap
            else:
                print("This attribute doesn't exist or is not implemented yet")
                break

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
                    if args.track_color == "u_tot":
                        c=e[j][:,8]
                    elif args.track_color == "beta_z":
                        c=e[j][:,4]/np.sqrt(1+e[j][:,8]**2)
                    elif args.track_color == "beta_y":
                        c=e[j][:,3]/np.sqrt(1+e[j][:,8]**2)
                    # print(len(x))
                    # print(len(y))
                    print("laenge track proc tag %i part tag %i ist %i" %(i, j, len(z)))
                    
                    # if args.twodproj:
                    #     plot_2D_colourline(z,x,c, cmin, cmax)
                    # else:
                    #     plot_3D_colourline(z,y,x,c, cmin, cmax)
                    
                    if args.twodproj:
                        if args.track_color == "u_tot":
                            plot_2D_colourline(z,x,c, cmin, cmax)
                        elif args.track_color == "beta_z":
                            plot_2D_colourline_beta(z,x,c)
                        elif args.track_color == "beta_y":
                            plot_2D_colourline_beta(z,x,c)
                    else:
                        if args.track_color == "u_tot":
                            plot_3D_colourline(z,y,x,c, cmin, cmax)
                        elif args.track_color == "beta_z":
                            plot_3D_colourline_beta(z,y,x,c)
                        elif args.track_color == "beta_y":
                            plot_3D_colourline_beta(z,y,x,c)
                    #ax.plot(z, y, x, label='particle track')
    
            ''' Set colorbar ''' 
            norm = matplotlib.colors.Normalize(
            vmin=np.min(cmin),
            vmax=np.max(cmax))

            # choose a colormap
            c_m = chosencmap

            # create a ScalarMappable and initialize a data structure
            s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
            s_m.set_array([])
    
            cbar = plt.colorbar(s_m)
            if args.track_color == "u_tot":
                cbar.ax.set_ylabel(r'$|u|$')
            elif args.track_color == "beta_z":
                cbar.ax.set_ylabel(r'$\beta_z$')
            elif args.track_color == "beta_y":
                cbar.ax.set_ylabel(r'$\beta_y$')
            
            
            
            ax.set_xlim(150, 200)
            ax.set_ylim(4, 6)
            ax.grid()
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