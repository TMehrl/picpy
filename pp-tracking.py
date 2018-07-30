#!/usr/bin/env python3
import csv
import os
import glob
import sys
import numpy as np
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import pp_defs
from pp_h5dat import mkdirs_if_nexist
from pp_h5dat import H5Plot
from mpl_toolkits.mplot3d import Axes3D
from pp_h5dat import Grid3d

from matplotlib.colors import LinearSegmentedColormap
import matplotlib.lines as mlines

import scipy.optimize
import scipy.special



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
                          action="store_true",
                          help = "Plots a 2D projection of the particle" 
                          "trajectories instead of full 3D plot. \n"
                          "Use only if acquired data is based on a single slice! (default: %(default)s).")

    parser.add_argument(  "-t", "--track_color",
                          dest = "track_color",
                          default = "u_tot",
                          help = "Select the attribute which is displayed as the color of the track \n"
                                 "Default is the total momentum u_tot, other available options are: \n"
                                 "beta_z")
    parser.add_argument(  "--h5",
                          dest = "h5plot",
                          action="store_true",
                          default=True,
                          help = "Save plot as hdf5 file (Default: %(default)s).")

    parser.add_argument(  "--dens_data",
                          dest = "dens_data",
                          default = "./DATA/density_ionized_electrons_plasma_H_000000.0.h5",
                          help = "Path to the ionized electron density data over which  \n"
                                 "the tracks should be laid. \n")

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
    # ax = fig.add_subplot(111)
    idx =np.abs(By_g3d2.get_zeta_arr() - zeta_pos).argmin()
    
    return By_g3d2.get_x_arr(2), Ex[:,idx]

    # savepath = './plots/g3d-line/'
    # mkdirs_if_nexist(savepath)
    # 
    # 
    # savename = 'Ex_analytic'
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # 
    # ax.plot(np.append(-r_array[::-1], r_array), np.append(-(E[::-1]) ,E[:]))
    # x, Ex =plot_hipace_Ex(-0.2) # input = zeta pos 
    # ax.plot(x, Ex)
    # 
    # h5lp = H5Plot()
    # h5lp.inherit_matplotlib_line_plots(ax)
    # h5lp.write(savepath + '/' + savename + '.h5')


def main():
    
    modnum = 6
    
    basic_cols=['#75b765', '#808080', '#ffd700']
    basic_cols=['#2020ff' ,'#808080', '#ff2020']
    basic_cols=['#0000ff' , '#00ffff','#808080', '#ffff00', '#ff0000']
    my_cmap=LinearSegmentedColormap.from_list('mycmap', basic_cols)
    


    parser = binSlab_parser() 
    args = parser.parse_args()

    
    dens_path = args.dens_data
    
    ionized_density_g3d1 = Grid3d(dens_path)
    ionized_density = np.transpose(ionized_density_g3d1.read(x2=0.0))


    # NHC=2
    proc_suffix_str = '_proc_'
    ppart_track_str = 'ppart_track'
    bin_fending = '.bin'
    
    
    zeta_min = -12

    zeta_max = 0.5
    zeta_gridpoints = 1200
    # zeta_min = -8
    # zeta_max = 4
    # zeta_gridpoints = 300
    zeta_array = np.arange(zeta_min, zeta_max, abs(zeta_max - zeta_min)/zeta_gridpoints)

    # print(datan)


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
        
        starting_positions = []
        
        for k in range(len(w)):
        
        
            fig = plt.figure()
            if args.twodproj:
                ax = fig.add_subplot(111)
            else:
                ax = fig.add_subplot(111, projection='3d')
            ''' get min and max value for universal colorbar later '''
            
            plt.pcolormesh(ionized_density_g3d1.get_zeta_arr(), ionized_density_g3d1.get_x_arr(2), ionized_density, cmap=cm.Blues) #
            c_m = cm.Blues
            s_m = matplotlib.cm.ScalarMappable(cmap=c_m)
            s_m.set_array([])
            cbar1 = plt.colorbar(s_m)
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
            # ##################### EXCLUDE HIPACE TRACKS FOR THE MOMENT
            ''' Splitting the array into each particle trajectory by proc tag'''
            d= np.split(w[k], np.where(np.diff(w[k][:, 7]) != 0)[0]+1) 
        
            for i in range(len(d)):
        
                ''' Splitting the array into each particle trajectory by part tag'''
                e = np.split(d[i], np.where(np.diff(d[i][:, 6]) != 0)[0]+1) 
                #print(np.shape(e))
                #print('particle tags %i (always look at 1 and ignore 0) %i'% (i, e[16][1, 6]))
                for j in range(int(np.floor(len(e)/modnum))):
                    x=e[modnum*j][:,0]
                    y=e[modnum*j][:,1]
                    z=e[modnum*j][:,5] # IN CASE OF UNIONIZED PLASMA USE THIS TERM
                    # starting_positions.append(x[ zeta_gridpoints -1]) #799])#
                    z=zeta_array[zeta_gridpoints-len(y):] # IN CASE OF preionized PLASMA USE THIS TERM
                    z=zeta_array[:len(y)]
                    if args.track_color == "u_tot":
                        c=e[modnum*j][:,8]
                    elif args.track_color == "beta_z":
                        c=e[modnum*j][:,4]/np.sqrt(1+e[modnum*j][:,8]**2)
                    elif args.track_color == "beta_y":
                        c=e[modnum*j][:,2]/np.sqrt(1+e[modnum*j][:,8]**2)
                    # print(len(x))
                    # print(len(y))
                    print("laenge track proc tag %i part tag %i ist %i" %(i, modnum*j, len(z)))
        
                    # if args.twodproj:
                    #     plot_2D_colourline(z,x,c, cmin, cmax)
                    # else:
                    #     plot_3D_colourline(z,y,x,c, cmin, cmax)
        
                    # ################### TAKEN OUT TO FASTEN EVERYTHING FOR THE 1/R analysis!
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
                    ###### ## ##ax.plot(z, y, x, label='particle track')
        
            # if args.twodproj:
            # 
            # 
            # 
            #     input_r_array = np.array([ 0, 0.0941176488995552, 0.1882352977991104, 0.2823529541492462, 0.3764705955982208, 0.47058823704719543, 0.5647059082984924, 0.658823549747467, 0.7529411911964417, 0.8470588326454163, 0.9411764740943909, 1.0352941751480103, 1.1294118165969849, 1.2235294580459595, 1.317647099494934, 1.4117647409439087, 1.5058823823928833, 1.600000023841858, 1.6941176652908325, 1.7882353067398071, 1.8823529481887817, 1.9764705896377563])
            #     input_r_array = np.array([0.02346041053533554, 0.04692082107067108, 0.07038123160600662, 0.09384164214134216, 0.1173020526766777, 0.14076246321201324, 0.16422288119792938, 0.18768328428268433, 0.21114370226860046, 0.2346041053533554, 0.25806450843811035, 0.2815249264240265])
            #     #input_r_array = np.array([0.3128054738044739, 0.703812301158905, 1.094819188117981, 1.485826015472412, 1.8768328428268433, 2.2678396701812744, 2.658846616744995, 3.0498533248901367, 3.4408602714538574, 3.831866979598999, 4.222873687744141, 4.613880634307861]) #rp 5 mod 1
            #     input_r_array = np.array([ 0.03519061580300331, 0.07038123160600662, 0.10557185113430023, 0.14076246321201324, 0.17595307528972626, 0.21114370226860046, 0.24633431434631348]) #rp 0.3 mod 3
            #     #input_r_array = np.array([0.033,0.034,0.03516, 0.03517,0.03519061580300331]) #,0.0352])
            #     # fig = plt.figure()
            #     # if args.twodproj:
            #     #     ax = fig.add_subplot(111)
            #     # else:
            #     #     ax = fig.add_subplot(111, projection='3d')
            #     # 
            #     # for i in range(int(np.floor(len(datan_zeta)/modnum))):
            #     #     ax.plot(datan_zeta[modnum*i,:], datan[modnum*i,:], color = '#00cc00', linestyle = '--') #'#551a8b'
            #     # # 
            #     # # 
            #     modnum = 1
            #     #for i in range(int(np.floor(len(input_r_array)/modnum))):
            #     for i in range(len(input_r_array)):
            #         zeta_array3, r_matrix2 = calc_analytical_solution(input_r_array[i], 5, 1, 2, 12, 0.3) #36, 0.3)#8,2)  ###(start_r_array, nb, ni, lbunch, zeta_end, rbunch):
            #         ax.plot(zeta_array3[np.where(r_matrix2 >0)], r_matrix2[np.where(r_matrix2 >0)], color = 'black',linestyle = '-.' ) ##00cc00
            #     # zeta_array2, r_matrix = calc_analytical_solution( 1.999999999, 1.2, 1, 2, 8,2) #1.8823529481887817 1.9764705896377563
            #     # ax.plot(zeta_array2, r_matrix, 'r' )  
            # else:
            #     ax = plt.gca()
            #     # ax.plot(data_zeta, np.zeros(len(data_zeta)), data2)
        
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
        
        
        
            #ax.set_xlim(-8, 0)
            # # ax.set_xlim(0, 300)
            #ax.set_ylim(-1/2, 6)
            # ax.grid()
            if not args.twodproj:
                ax.set_zlabel(' x ')
                #ax.set_zlim(-6, 6)
            ax.set_xlabel(r'$\zeta$')
            ax.set_ylabel(' y ')
        
        
            print(starting_positions)
            
            #plt.show()
            # fname, fext = os.path.splitext(args.paths[0])


            save_path_name = 'plots/g3d-slice/ionized_electron_density_tracked.pdf'
            print('bis kury vor dem safe')
            fig.savefig(save_path_name, format='png')
            print('bis nach dem save')
        #    plt.close(fig)
            if args.verbose: 
                sys.stdout.write('Saved: %s\n' % save_path_name)
                sys.stdout.flush()
            plt.close(fig)
            print('bis nach dem close fig')
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