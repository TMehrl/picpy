#!/usr/bin/env python3
# This script may be executed like this:
# nohup ./raw-slice-series-plotting.py <DATA>/ 1> rss.out 2> rss.err &

import os
import sys
import argparse

import math
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib.colors import LogNorm
from matplotlib.ticker import FormatStrFormatter
import pp_defs
import h5py
from pp_h5dat import SliceMoms
from pp_h5dat import H5Plot
from pp_h5dat import mkdirs_if_nexist
import pp_raw_ana


# Parse defaults/definitions
class parsedefs:
    class file_format:
        png = 'png'
        eps = 'eps'
        pdf = 'pdf'
    class zax:
        zeta = 'zeta'
        z = 'z'
        xi = 'xi'
    class save:
        prefix = 'rss_name'
        path = './plots'

def ps_parseopts():

    desc="""This is the picpy postprocessing tool."""

    parser = argparse.ArgumentParser(description=desc)
    #parser = OptionParser(usage=usg, description=desc)
    parser.add_argument(  'path',
                          metavar = 'PATH',
                          help = 'Path to raw files.')    
    parser.add_argument(  "-v", "--verbose",
                        action="store_true",
                        dest="verbose",
                        default=True,
                        help = "Print info (default: %(default)s).")
    parser.add_argument(  "-q", "--quiet",
                        action="store_false",
                        dest="verbose",
                        help = "Don't print info.")
    parser.add_argument(  "-s", "--save-path",
                          action="store",
                          dest="savepath",
                          metavar="PATH",
                          default=parsedefs.save.path + '/raw-slice-series',
                          help = """Path to which generated files will be saved.
                              (default: %(default)s)""")
    parser.add_argument(  "-n", "--name-prefix",
                        dest="save_prefix",
                        metavar="NAME",
                        default=parsedefs.save.prefix,
                        help = """Define customized prefix of output filename.""")
    parser.add_argument(  "-c", "--code",
                        action='store',
                        dest="piccode",
                        metavar="CODE",
                        choices = [pp_defs.code.hipace, pp_defs.code.osiris,],
                        default = pp_defs.code.hipace,
                        help= "PIC code which was used to generate files (default: %(default)s).")
    parser.add_argument(  "-d", "--dim",
                        action='store',
                        dest="dimensionality",
                        metavar="DIM",
                        choices=[1, 2, 3,],
                        default=3,
                        help= """Dimensionality of PIC simulation
                              (default: %(default)s).""")
    parser.add_argument(  "-N", "--number-of-files",
                        action='store',
                        dest="Nfiles",
                        metavar="NFILES",
                        default=0,
                        help= """Number of files to analyze.""")
    parser.add_argument(  "--h5",
                          dest = "h5plot",
                          action="store_true",
                          default=True,
                          help = "Save plot as hdf5 file (Default: %(default)s).")    
    return parser


def plot_save_slice_rms(slm, savepath):

    x = slm.zeta_array
    y = slm.time_array

    sigma_x = np.sqrt( np.absolute( slm.avgx2sq ) )
    fig_sx = plt.figure()
    cax = plt.pcolormesh( x,
                          y,
                          sigma_x,
                          cmap=cm.Blues,
                          vmin=0, vmax=np.amax(abs(sigma_x)) )
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$\omega_p t$', fontsize=14)    
    cbar = fig_sx.colorbar( cax )
    cbar.ax.set_ylabel(r'$k_p \sigma_x$', fontsize=14)
    fig_sx.savefig( savepath + '/sigma_x.png',
                  format='png',
                  dpi=600)


    sigma_px = np.sqrt( np.absolute( slm.avgp2sq ) )


    fig_spx = plt.figure()
    cax = plt.pcolormesh( x,
                          y,
                          sigma_px,
                          cmap=cm.YlGn,
                          vmin=0, vmax=np.amax(sigma_px) )
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$\omega_p t$', fontsize=14)    
    cbar = fig_spx.colorbar( cax )
    cbar.ax.set_ylabel(r'$k_p \sigma_{p_x}$', fontsize=14)
    fig_spx.savefig( savepath + '/sigma_px.png',
                  format='png',
                  dpi=600)

    emittance = np.sqrt( np.multiply(slm.avgx2sq, slm.avgp2sq) 
                         - np.power(slm.avgx2p2,2) )
    fig_e = plt.figure()
    cax = plt.pcolormesh( x,
                          y,
                          emittance,
                          cmap=cm.Reds,
                          vmin=np.amin(emittance), vmax=np.amax(emittance) )
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$\omega_p t$', fontsize=14)    
    cbar = fig_e.colorbar( cax )
    cbar.ax.set_ylabel(r'$k_p \epsilon_x$', fontsize=14)
    fig_e.savefig( savepath + '/slice_emittance_x.png',
                  format='png',
                  dpi=600)


def plot_save_proj_rms(slm, savepath, h5plot=True):
    
    t = slm.time_array
    tot_charge = np.sum(slm.charge, axis=1)
    xsq = np.divide(np.sum(np.multiply(slm.avgx2sq, slm.charge), axis=1),tot_charge)
    psq = np.divide(np.sum(np.multiply(slm.avgp2sq, slm.charge), axis=1),tot_charge)
    xp = np.divide(np.sum(np.multiply(slm.avgx2p2, slm.charge), axis=1),tot_charge)


    fig_sx = plt.figure()    
    plt.plot(t, np.sqrt(xsq))
    ax = plt.gca()
    ax.set_xlabel(r'$\omega_p t$', fontsize=14) 
    ax.set_ylabel(r'$k_p\sigma_x$', fontsize=14)
    if not (-3.0 < math.log(np.max(abs(np.sqrt(xsq))),10) < 3.0):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15) 
    fig_sx.savefig( savepath + '/sigma_x_proj.eps',
                  format='eps')
    if h5plot: 
        h5lp = H5Plot()
        h5lp.inherit_matplotlib_line_plots(ax)
        h5lp.write(savepath + '/sigma_x_proj.h5')
    plt.close(fig_sx)


    fig_sp = plt.figure()    
    plt.plot(t, np.sqrt(psq))
    ax = plt.gca()
    ax.set_xlabel(r'$\omega_p t$', fontsize=14) 
    ax.set_ylabel(r'$\sigma_{p_x}/m_e c$', fontsize=14)
    if not (-3.0 < math.log(np.max(abs(np.sqrt(psq))),10) < 3.0):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15)  
    fig_sx.savefig( savepath + '/sigma_px_proj.eps',
                  format='eps')
    if h5plot: 
        h5lp = H5Plot()
        h5lp.inherit_matplotlib_line_plots(ax)
        h5lp.write(savepath + '/sigma_px_proj.h5')
    plt.close(fig_sp)


    fig_xp = plt.figure()    
    plt.plot(t, xp)
    ax = plt.gca()
    ax.set_xlabel(r'$\omega_p t$', fontsize=14) 
    ax.set_ylabel(r'$ k_p x\,p_x/m_e c$', fontsize=14)
    if not (-3.0 < math.log(np.max(abs(xp)),10) < 3.0):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15)
    fig_xp.savefig( savepath + '/xpx_proj.eps',
                  format='eps')
    if h5plot: 
        h5lp = H5Plot()
        h5lp.inherit_matplotlib_line_plots(ax)
        h5lp.write(savepath + '/xpx_proj.h5')
    plt.close(fig_xp)


    emittance = np.sqrt(xsq*psq-np.power(xp,2))
    fig_e = plt.figure()    
    plt.plot(t, emittance)
    ax = plt.gca()
    ax.set_xlabel(r'$\omega_p t$', fontsize=14) 
    ax.set_ylabel(r'$k_p \epsilon_x$', fontsize=14)
    if not (-3.0 < math.log(np.max(abs(emittance)),10) < 3.0):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15)  
    fig_e.savefig( savepath + '/emittance_proj.eps',
                  format='eps')
    if h5plot: 
        h5lp = H5Plot()
        h5lp.inherit_matplotlib_line_plots(ax)
        h5lp.write(savepath + '/emittance_proj.h5')
    plt.close(fig_e)

def plot_save_slice_centroids(slm, savepath):

    Xb0 = np.ones(slm.avgx2[0,:].shape)
    Yb0 = np.ones(slm.avgx3[0,:].shape)
    for i in range(0,len(slm.zeta_array)):
        if (slm.avgx2[0,i] != 0.0):
            Xb0[i] = slm.avgx2[0,i]
            Yb0[i] = slm.avgx3[0,i]


    Xb_norm = np.zeros( slm.avgx2.shape )
    Yb_norm = np.zeros( slm.avgx3.shape )
#    zeta_hseed = np.min(slm.zeta_array)
#    idx_hseed = (np.abs(slm.zeta_array-zeta_hseed)).argmin()

    for i in range(0,len(slm.zeta_array)):
#        if (slm.zeta_array[i] <= zeta_hseed):
        Xb_norm[:,i] = np.absolute( slm.avgx2[:,i]/Xb0[i] )
        Yb_norm[:,i] = np.absolute( slm.avgx3[:,i]/Yb0[i] )

    figXb = plt.figure()
    cax = plt.pcolormesh(   slm.zeta_array,
                            slm.time_array,
                            slm.avgx2,
                            cmap=cm.PuOr,
                            vmin=-np.amax(abs(slm.avgx2)), 
                            vmax=np.amax(abs(slm.avgx2)))
    cbar = figXb.colorbar(cax)   
    cbar.ax.set_ylabel('$k_p X_b$', fontsize=14)
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$\omega_p t$', fontsize=14)     
    figXb.savefig( savepath + '/Xb_raw.png',
                  format='png',
                  dpi=600)

    figYb = plt.figure()
    cax = plt.pcolormesh(   slm.zeta_array,
                            slm.time_array,
                            slm.avgx3,
                            cmap=cm.PuOr,
                            vmin=-np.amax(abs(slm.avgx3)), 
                            vmax=np.amax(abs(slm.avgx3)))
    cbar = figYb.colorbar(cax)
    cbar.ax.set_ylabel('$k_p Y_b$', fontsize=14)
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$\omega_p t$', fontsize=14)     
    figYb.savefig( savepath + '/Yb_raw.png',
                  format='png',
                  dpi=600)

    figXbnorm = plt.figure()
    cax = plt.pcolormesh(   slm.zeta_array,
                            slm.time_array,
                            Xb_norm,
                            cmap=cm.Blues,
                            vmin=0, 
                            vmax=np.amax(Xb_norm))
    cbar = figXbnorm.colorbar(cax)
    cbar.ax.set_ylabel('$|X_b/X_{b,0}|$', fontsize=14)
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$\omega_p t$', fontsize=14)     
    figXbnorm.savefig( savepath + '/Xb.png',
                  format='png',
                  dpi=600)

    figYbnorm = plt.figure()
    cax = plt.pcolormesh(   slm.zeta_array,
                            slm.time_array,
                            Yb_norm,
                            cmap=cm.Blues,
                            vmin=0, 
                            vmax=np.amax(Yb_norm))
    cbar = figYbnorm.colorbar(cax)
    cbar.ax.set_ylabel('$|Y_b/Y_{b,0}|$', fontsize=14)
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$\omega_p t$', fontsize=14)       
    figYbnorm.savefig( savepath + '/Yb.png',
                  format='png',
                  dpi=600)

    figXb0 = plt.figure()
    plt.plot(slm.zeta_array, Xb0)
    ax = plt.gca()
    ymin, ymax = ax.get_ylim()
    if ymin > 0 and ymax > 0:
        plt.ylim(0, ymax*1.2)
    elif ymin < 0 and ymax < 0:
        plt.ylim(ymin*1.2, 0)        
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$k_p X_{b,0}$', fontsize=14)       
    figXb0.savefig( savepath + '/Xb0.eps',
                  format='eps')


    figYb0 = plt.figure()
    plt.plot(slm.zeta_array, Yb0)
    ax = plt.gca()
    ymin, ymax = ax.get_ylim()
    if ymin > 0 and ymax > 0:
        plt.ylim(0, ymax*1.2)
    elif ymin < 0 and ymax < 0:
        plt.ylim(ymin*1.2, 0)    
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$k_p Y_{b,0}$', fontsize=14)     
    figYb0.savefig( savepath + '/Yb0.eps',
                  format='eps')


    figXbtail = plt.figure()
    plt.plot(slm.time_array, slm.avgx2[:,0])
    ax = plt.gca()
    ax.set_xlabel(r'$\omega_p t$', fontsize=14)
    ax.set_ylabel(r'$X_{b,\mathrm{tail}}$', fontsize=14)
    if not (-3.0 < math.log(np.max(abs(slm.avgx2[:,0])),10) < 3.0):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15)
    figXbtail.savefig( savepath + '/Xb_tail.eps',
                  format='eps')

    figYbtail = plt.figure()
    plt.plot(slm.time_array, slm.avgx3[:,0])
    ax = plt.gca()
    ax.set_xlabel(r'$\omega_p t$', fontsize=14)
    ax.set_ylabel(r'$Y_{b,\mathrm{tail}}$', fontsize=14) 
    figYbtail.savefig( savepath + '/Yb_tail.eps',
                  format='eps')


def plot_save_slice_ene(slm, savepath):
    gamma = np.sqrt( 1 + np.power(slm.avgp1,2) 
                       + np.power(slm.avgp2,2)
                       + np.power(slm.avgp3,2) )

    fig7 = plt.figure()
    cax = plt.pcolormesh( slm.zeta_array,
                          slm.time_array,
                          gamma,
                          cmap=cm.GnBu,
                          vmin=0, vmax=np.amax(gamma))
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$\omega_p t$', fontsize=14)    
    cbar = fig7.colorbar( cax )
    cbar.ax.set_ylabel(r'$\gamma$', fontsize=14)
    fig7.savefig( savepath + '/gamma.png',
                  format='png',
                  dpi=600)

def plot_curr_profile(slm, savepath):
    fig = plt.figure()
    dzeta = abs(slm.zeta_array[1] - slm.zeta_array[0]);
    curr = slm.charge[0,:] / dzeta

    print('Q = %0.3e' % (np.sum(curr) * dzeta ))

    # I_A = 4 * pi * epsilon_0 * m * c^3 / e
    # [curr] = epsilon_0 * m * c^3 / e
    # curr * 4 * pi = I_b/I_A
    plt.plot(slm.zeta_array, curr/(4 * math.pi))
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$I_{b,0}/I_A$', fontsize=14) 
    fig.savefig( savepath + '/Ib0.eps',
                  format='eps')


def main():

    parser = ps_parseopts()

    args = parser.parse_args()

    file = args.path

    slm = SliceMoms(file)

    mkdirs_if_nexist(args.savepath)

    plot_curr_profile(slm, args.savepath)

    plot_save_slice_rms(slm, args.savepath)

    plot_save_slice_ene(slm, args.savepath)

    plot_save_proj_rms(slm, args.savepath, args.h5plot)

    plot_save_slice_centroids(slm, args.savepath)


if __name__ == "__main__":
    main()
