#!/usr/bin/env python3
import csv
import os
import sys
import math
import scipy.optimize
import scipy.special
import numpy as np
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.lines as mlines
from matplotlib.ticker import MaxNLocator


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
    parser.add_argument(  "--save-name",
                          action="store",
                          dest="savename",
                          metavar="SNAME",
                          default='/ionized_electron_density_',
                          help = """Name of the file + eventually 'tracked' + timestamp""")
    parser.add_argument(  "--nohalo",
                          dest = "nohalo",
                          action="store_true",
                          default=True,
                          help = "If data is w/o halo (default: %(default)s).")
    parser.add_argument(  "-t", "--track_color",
                          dest = "track_color",
                          default = "none",
                          help = "Select the attribute which is displayed as the color of the track \n"
                                 "Default is the total momentum u_tot, other available options are: \n"
                                 "beta_z")
    parser.add_argument(  "--dens_data",
                          dest = "dens_data",
                          nargs = '*',
                          default = "./DATA/density_ionized_electrons_plasma_H_000000.0.h5",
                          help = "Path to the ionized electron density data over which  \n"
                                 "the tracks should be laid. \n")
    parser.add_argument(  '--clim',
                          help='Colorbar axis limits',
                          action='store',
                          dest="clim",
                          metavar=('CBMIN', 'CBMAX'),
                          nargs=2,
                          type=float,
                          default=None)
    parser.add_argument(  "--beam_data",
                        dest = "beam_data",
                        nargs = '*',
                        default = False,
                        help = "Path to the beam density data over which  \n"
                               "the tracks should be laid. \n")
    parser.add_argument(  '--cblim',
                        help='Colorbar axis limits for beam density colorbar',
                        action='store',
                        dest="cblim",
                        metavar=('CBMIN', 'CBMAX'),
                        nargs=2,
                        type=float,
                        default=None)
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
                          default=None)
    parser.add_argument(  "--manual",
                          dest = "manual",
                          action="store_true",
                          default=False,
                          help = "Enables manual multiplots (timestamps don't have to agree) (Default: %(default)s).")
    parser.add_argument(  "--old_timestamp",
                          dest = "old_timestamp",
                          action="store_true",
                          default=False,
                          help = "Enables tracking of TRACKDATA with old timestamp (Default: %(default)s).")
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
    parser.add_argument(  '--lw',
                        help='linewidth of the particle tracks',
                        action='store',
                        dest="linewidth",
                        type=float,
                        default=0.3)
    parser.add_argument(  "--tracksoff",
                          dest = "tracksoff",
                          action="store_true",
                          default=False,
                          help = "Dismiss particle tracks (Default: %(default)s).")
    parser.add_argument(  "--ptype",
                          default="contourf",
                          dest="ptype",
                          choices=[ "pcolormesh", "contourf"],
                          help= "Plot color type (default: %(default)s).")
    parser.add_argument(  "--latexfont",
                          dest = "latexfont",
                          action="store_true",
                          default=False,
                          help = "Use LaTeX font (Default: %(default)s).")

    return parser



def plot_2D_colourline(x,z,c, minc, maxc, lw):

    c = cm.jet((c-minc)/(maxc-minc))
    ax = plt.gca()
    for i in np.arange(len(x)-1):
        ax.plot([x[i],x[i+1]], [z[i],z[i+1]], c=c[i], linewidth=lw)

    return


def plot_2D_colourline_beta(x,z,c, lw):
    basic_cols=['#2020ff' ,'#808080', '#ff2020']
    basic_cols=['#0000ff' , '#00ffff','#808080', '#ffff00', '#ff0000']
    my_cmap = LinearSegmentedColormap.from_list('mycmap', basic_cols)
    c = my_cmap((c+1)/(2))
    ax = plt.gca()
    for i in np.arange(len(x)-1):
        ax.plot([x[i],x[i+1]], [z[i],z[i+1]], c=c[i], linewidth=lw)

    return

def round_figures(x, n):
    """Returns x rounded to n significant figures."""
    return round(x, int(n - math.ceil(math.log10(abs(x)))))

def main():



    basic_cols=['#75b765', '#808080', '#ffd700']
    basic_cols=['#2020ff' ,'#808080', '#ff2020']
    basic_cols=['#0000ff' , '#00ffff','#808080', '#ffff00', '#ff0000']
    my_cmap=LinearSegmentedColormap.from_list('mycmap', basic_cols)

    parser = binSlab_parser()
    args = parser.parse_args()

    if args.latexfont:
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        
    modnum = args.modlines
    density_path = args.dens_data


    if not args.tracksoff and not args.track_range:
        print('ERROR: tracks not turned off, but no tracking range given. Provide tracking range with --track_range to get correct colormap.')

    # NHC=2
    proc_suffix_str = '_proc_'
    ppart_track_str = 'ppart_track'
    bin_fending = '.bin'



    for density_files in args.dens_data:

        ionized_density_g3d1 = Grid3d(density_files)
        ionized_density = np.transpose(ionized_density_g3d1.read(x2=0.0))
        print(density_files)
        timestamp = density_files.split("_")[-1].split(".h5")[0]
        if args.old_timestamp:
            timestamp = timestamp.split(".")[0]

        if args.beam_data:
            for beam_density_files in args.beam_data:
                if timestamp in beam_density_files:
                    beam_density_g3d1 = Grid3d(beam_density_files)
                    beam_density = np.transpose(beam_density_g3d1.read(x2=0.0))

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
                        if (timestamp in files) or args.manual:
                            print('Reading: %s/%s' % (dirname, files))
                            filepath = dirname + slash + files
                            array = np.append(array,np.fromfile(filepath,dtype=np.float32))
                else:
                    if (timestamp in files) or args.manual:
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

                # if not args.tracksoff and args.dens_data and args.beam_data:
                #     fig = plt.figure(figsize=(9,5))
                # else:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ''' get min and max value for universal colorbar later '''
                max_density = np.max(np.abs(ionized_density))
                min_density = 0
                if args.clim:
                    levels = MaxNLocator(nbins=512).tick_values(args.clim[0], args.clim[1])
                    max_level = max_density
                    vmin = args.clim[0]
                    vmax = args.clim[1]
                    #selecting correct extend method
                    if min_density < args.clim[0] and max_density > args.clim[1]:
                        extend = 'both'
                    elif min_density < args.clim[0] and max_density <= args.clim[1]:
                        extend = 'min'
                    elif min_density >= args.clim[0] and max_density > args.clim[1]:
                        extend = 'max'
                    elif min_density >= args.clim[0] and max_density <= args.clim[1]:
                        extend = 'neither'
                    else:
                        print('Error: unexpected case, couldn\'t extend in the correct way!')
                        extend = 'neither'

                else:
                    levels = MaxNLocator(nbins=512).tick_values(0, max_density)
                    max_level = max_density
                    vmin = min_density
                    vmax = max_density
                    extend = 'neither'
                if args.ptype == 'pcolormesh':
                    plot1 = plt.pcolormesh(ionized_density_g3d1.get_zeta_arr(), ionized_density_g3d1.get_x_arr(2), np.abs(ionized_density), cmap=cm.PuBu, alpha=1) #
                elif args.ptype == 'contourf':
                    plot1 = plt.contourf(ionized_density_g3d1.get_zeta_arr(), ionized_density_g3d1.get_x_arr(2), np.abs(ionized_density), cmap=cm.PuBu, levels=levels, vmin=vmin, vmax=vmax, extend=extend)
                else:
                    print('This type is not implemented yet')

                cbarmap = plt.cm.ScalarMappable(cmap=cm.PuBu)
                cbarmap.set_array(np.abs(ionized_density))
                if args.clim:
                    cbarmap.set_clim(args.clim[0], args.clim[1])
                    cbar1= plt.colorbar(cbarmap, boundaries=np.arange(args.clim[0],args.clim[1]+0.0002,0.0001), extend=extend, fraction=0.046, pad=0.1)
                    ticks = MaxNLocator(5).tick_values(vmin, vmax)
                    cbar1.set_ticks ( ticks )
                    plot1.set_clim([args.clim[0], args.clim[1]])

                else:
                    cbarmap.set_clim([0, max_density])
                    cbar1= plt.colorbar(cbarmap, boundaries=np.arange(0,max_density+0.0002,0.0001), extend=extend, fraction=0.046, pad=0.1)
                    ticks = MaxNLocator(5).tick_values(vmin, vmax)
                    cbar1.set_ticks ( ticks )
                    plot1.set_clim([0, max_density])

                if args.xlim:
                    plt.xlim(args.xlim[0], args.xlim[1])
                if args.ylim:
                    plt.ylim(args.ylim[0], args.ylim[1])


                cbar1.ax.set_title(r'$n_p/n_0$')

                if args.beam_data:
                    max_density = np.max(np.abs(beam_density))

                    cmap = plt.cm.Reds

                    # Get the colormap colors
                    my_cmap = cmap(np.arange(cmap.N))

                    # Set alpha
                    my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
                    # Create new colormap
                    my_cmap = ListedColormap(my_cmap)

                    if args.cblim:
                        levels = MaxNLocator(nbins=512).tick_values(args.cblim[0], args.cblim[1])
                        max_level = max_density
                        vmin = args.cblim[0]
                        vmax = args.cblim[1]
                        #selecting correct extend method
                        if min_density < args.cblim[0] and max_density > args.cblim[1]:
                            extend = 'both'
                        elif min_density < args.cblim[0] and max_density <= args.cblim[1]:
                            extend = 'min'
                        elif min_density >= args.cblim[0] and max_density > args.cblim[1]:
                            extend = 'max'
                        elif min_density >= args.cblim[0] and max_density <= args.cblim[1]:
                            extend = 'neither'
                        else:
                            print('Error: unexpected case, couldn\'t extend in the correct way!')
                            extend = 'neither'
                    else:
                        levels = MaxNLocator(nbins=512).tick_values(0, max_density)
                        max_level = max_density
                        vmin = 0
                        vmax = max_density
                        extend = 'neither'
                    if args.ptype == 'pcolormesh':
                        plt.pcolormesh(beam_density_g3d1.get_zeta_arr(), beam_density_g3d1.get_x_arr(2), np.abs(beam_density), cmap=my_cmap, vmin=vmin, vmax=vmax)
                    elif args.ptype == 'contourf':
                        plt.contourf(beam_density_g3d1.get_zeta_arr(), beam_density_g3d1.get_x_arr(2), np.abs(beam_density), cmap=my_cmap, levels=levels, vmin=vmin, vmax=vmax, extend=extend)
                    else:
                        print('This type is not implemented yet')


                    cbarmap = plt.cm.ScalarMappable(cmap=my_cmap)
                    cbarmap.set_array(np.abs(beam_density))
                    if args.cblim:
                        # Note on colorbar: boundaries have to be set manually, because otherwise there will be ugly stripes
                        # afterwards the ticks have to set manually as well, set them at the correct place
                        cbarmap.set_clim(args.cblim[0], args.cblim[1])
                        cbar2 = plt.colorbar(cbarmap, boundaries=np.arange(args.cblim[0],args.cblim[1]+0.0002,0.0001), extend=extend, fraction=0.046, pad=0.1 )
                        ticks = MaxNLocator(5).tick_values(vmin, vmax)
                        cbar2.set_ticks ( ticks )
                        cbar2.set_clim([args.cblim[0], args.cblim[1]])
                    else:
                        #cbarmap.set_clim(0, max_density)
                        cbar2= plt.colorbar(cbarmap, boundaries=np.arange(0,max_density+0.0001,0.0001), fraction=0.046, pad=0.1 )
                        ticks = MaxNLocator().tick_values(vmin, vmax)
                        cbar2.set_ticks ( ticks )
                        cbar2.set_clim([0, max_density])

                    cbar2.ax.set_title(r'$n_b/n_0$')


                if not args.tracksoff:

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
                                plot_2D_colourline(z,x,c, cmin, cmax, args.linewidth)
                            elif args.track_color == "beta_z":
                                plot_2D_colourline_beta(z,x,c, args.linewidth)
                            elif args.track_color == "beta_y":
                                plot_2D_colourline_beta(z,x,c, args.linewidth)
                            elif args.track_color == 'none':
                                if len(x)>0:
                                    index = np.argmin(abs(start_segments - x[-1]))

                                    # print('index : %i' %index)
                                    # print('x0 ist %f' %x[-1])
                                    ax.plot(z, x, color=colors2[index], linewidth = 0.3)



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

                        cbar = plt.colorbar(s_m, fraction=0.052, pad=0.06)

                    else:
                        # norm = matplotlib.colors.Normalize( vmin= 0, vmax=6)

                        # choose a colormap
                        c_m = cm.jet

                        # create a ScalarMappable and initialize a data structure
                        s_m = matplotlib.cm.ScalarMappable(cmap=c_m)
                        s_m.set_array([args.track_range[0],args.track_range[1]])

                        cbar = plt.colorbar(s_m, fraction=0.052, pad=0.06)

                    if args.track_color == "u_tot":
                        cbar.ax.set_title(r'$|u|$')
                    elif args.track_color == "beta_z":
                        cbar.ax.set_title(r'$\beta_z$')
                    elif args.track_color == "beta_y":
                        cbar.ax.set_title(r'$\beta_y$')
                    else:
                        cbar.ax.set_title(r'$k_p\,x_0$')

                ax.set_xlabel(r'$k_p\,\zeta$', fontsize=14)
                ax.set_ylabel(r'$k_p\,x$', fontsize=14)


                savepath = 'plots/g3d-slice'
                mkdirs_if_nexist(savepath)
                if not args.tracksoff:
                    tracked = 'tracked_'
                else:
                    tracked = ''
                if args.savename:
                    savename = '/' + args.savename + '_'
                else:
                    savename = '/ionized_electron_density_'
                save_path_name = savepath + savename + tracked + timestamp + '.png'
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
