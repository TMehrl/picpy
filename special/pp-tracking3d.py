#!/usr/bin/env python3
import csv
import os
import sys
import numpy as np
import scipy.optimize
import scipy.special
import argparse
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.lines as mlines

mypath = os.path.dirname(os.path.abspath( __file__ ))
incpath = os.path.split(mypath)[0] + '/inc'
sys.path.append(incpath)
import pp_defs
from pp_h5dat import mkdirs_if_nexist
from pp_h5dat import H5Plot
from mpl_toolkits.mplot3d import Axes3D
from pp_h5dat import Grid3d


def binSlab_parser():

    desc = """This is the picpy postprocessing tool."""

    parser = argparse.ArgumentParser( description=desc)
    parser.add_argument(  'path',
                          metavar = 'PATH',
                          nargs = '*',
                          help = 'Path to tracking files. Important: In order to use multiple tracking files, \n' 
                          'DO NOT use * just leave the end open like tracking_data_proc_ (instead of tracking_data_proc_* )')
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
    parser.add_argument(  "-t", "--track_color",
                          dest = "track_color",
                          default = "none",
                          help = "Select the attribute which is displayed as the color of the track \n"
                                 "Default is the total momentum u_tot, other available options are: \n"
                                 "beta_z")
    parser.add_argument(  '--modlines',
                        help='number of lines modulo input',
                        action='store',
                        dest="modlines",
                        type=int,
                        default=1) 
    parser.add_argument(  '--track_range',
                          help='Give the length of the original tracking box',
                          action='store',
                          dest="track_range",
                          metavar=('trackmin', 'trackmax'),
                          nargs=2,
                          type=float,
                          default=[0,1])
    parser.add_argument(  "--xlim",
                          help='Customize x-axis limits',
                          action='store',
                          dest="xlim",
                          metavar=('xmin', 'xmax'),
                          type=float,
                          nargs=2,
                          default=None)  
    parser.add_argument(  "--ylim",
                          help='Customize y-axis limits',
                          action='store',
                          dest="ylim",
                          metavar=('ymin', 'ymax'),
                          type=float,
                          nargs=2,
                          default=None)  


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
        ax.plot([x[i],x[i+1]], [z[i],z[i+1]], c=c[i], linewidth=0.2)
    
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


def plot_hipace_Ex(zeta_pos):
    ExmBy_path = './DATA/field_ExmBy_000000.h5'
    By_path = './DATA/field_By_000000.h5'
    
    ExmBy_g3d2 = Grid3d(ExmBy_path)
    By_g3d2 = Grid3d(By_path)
    
    ExmBy = np.transpose(ExmBy_g3d2.read(x2=0.0))
    By = np.transpose(By_g3d2.read(x2=0.0))  

    Ex = ExmBy + By
    
    fig = plt.figure()

    idx =np.abs(By_g3d2.get_zeta_arr() - zeta_pos).argmin()
    
    return By_g3d2.get_x_arr(2), Ex[:,idx]




def main():
    
    
    
    basic_cols=['#75b765', '#808080', '#ffd700']
    basic_cols=['#2020ff' ,'#808080', '#ff2020']
    basic_cols=['#0000ff' , '#00ffff','#808080', '#ffff00', '#ff0000']
    my_cmap=LinearSegmentedColormap.from_list('mycmap', basic_cols)
    
    parser = binSlab_parser() 
    args = parser.parse_args()


    modnum = args.modlines

    



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
            elif args.track_color == "none":
                print('no coloring option selected.')
            else:
                print("This attribute doesn't exist or is not implemented yet")
                break

            ''' Splitting the array into each particle trajectory by proc tag'''
            d= np.split(w[k], np.where(np.diff(w[k][:, 7]) != 0)[0]+1) 
        
            for i in range(len(d)):
        
                ''' Splitting the array into each particle trajectory by part tag'''
                e = np.split(d[i], np.where(np.diff(d[i][:, 6]) != 0)[0]+1) 

                
                number = int(np.floor(len(e)/modnum))
                cmap = plt.get_cmap('jet')
                colors = [cmap(i) for i in np.linspace(0, 1, number)]
                colors2 = [cmap(i) for i in np.linspace(0, 1, 10000)]
                start_segments = np.linspace(args.track_range[0],args.track_range[1],10000)

                
                for j in range(int(np.floor(len(e)/modnum))):
                    x=e[modnum*j][:,0]
                    y=e[modnum*j][:,1]
                    z=e[modnum*j][:,5] 
                    
                    if args.track_color == "u_tot":
                        c=e[modnum*j][:,8]
                    elif args.track_color == "beta_z":
                        c=e[modnum*j][:,4]/np.sqrt(1+e[modnum*j][:,8]**2)
                    elif args.track_color == "beta_y":
                        c=e[modnum*j][:,2]/np.sqrt(1+e[modnum*j][:,8]**2)
                    print("laenge track proc tag %i part tag %i ist %i" %(i, modnum*j, len(z)))

                    

                    if args.track_color == "u_tot":
                        plot_3D_colourline(z,x,y,c, cmin, cmax)
                    elif args.track_color == "beta_z":
                        plot_3D_colourline_beta(z,x,y,c)
                    elif args.track_color == "beta_y":
                        plot_3D_colourline_beta(z,x,y,c)
                    else: 
                        plot_3D_colourline(z,x,y,c, cmin, cmax)

        
            ''' Set colorbar ''' 
            if args.track_color != "none":
                norm = matplotlib.colors.Normalize(
                vmin=np.min(cmin),
                vmax=np.max(cmax))
            
                # choose a colormap
                c_m = chosencmap
            
                # create a ScalarMappable and initialize a data structure
                s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
                s_m.set_array([])
            
                cbar = plt.colorbar(s_m)
                
            else:
                # norm = matplotlib.colors.Normalize( vmin= 0, vmax=6)
            
                # choose a colormap
                c_m = cm.jet
            
                # create a ScalarMappable and initialize a data structure
                s_m = matplotlib.cm.ScalarMappable(cmap=c_m)
                s_m.set_array([args.track_range[0],args.track_range[1]])
            
                cbar = plt.colorbar(s_m)
                
            if args.track_color == "u_tot":
                cbar.ax.set_title(r'$|u|$', fontsize=18)
            elif args.track_color == "beta_z":
                cbar.ax.set_title(r'$\beta_z$', fontsize=18)
            elif args.track_color == "beta_y":
                cbar.ax.set_title(r'$\beta_y$', fontsize=18)
            else:
                cbar.ax.set_title(r'$k_p\,x_0$', fontsize=18)
        
        
                    #ax.set_xlim(-8, 0)
            # # ax.set_xlim(0, 300)
            #ax.set_ylim(-1/2, 6)
            if args.xlim:
                ax.set_ylim(args.xlim[0], args.xlim[1])
            if args.ylim:
                ax.set_zlim(args.ylim[0], args.ylim[1])

            ax.set_zlabel(r'$k_p\,y$', fontsize=18)
            ax.set_xlabel(r'$k_p\,\zeta$', fontsize=18)
            ax.set_ylabel(r'$k_p\,x$', fontsize=18)

            plt.show()
            
            savepath = 'plots/g3d-slice'
            mkdirs_if_nexist(savepath)
                    
            save_path_name = savepath + '/ionized_electron_density_tracked' + '.png'
            print('Saving figure...')
            fig.savefig(save_path_name, format='png', dpi=400)
            print('Writing file...')

            if args.verbose: 
                sys.stdout.write('Saved: %s\n' % save_path_name)
                sys.stdout.flush()
            plt.close(fig)
            print('Done!')


if __name__ == "__main__":
    main()  