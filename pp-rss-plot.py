#!/usr/bin/env python3
# This script may be executed like this:
# nohup ./raw-slice-series-plotting.py <DATA>/ 1> rss.out 2> rss.err &

import os
import sys
import argparse
from optparse import OptionParser
#from optparse import OptionGroup
import math
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib.colors import LogNorm
import pp_defs
from pp_h5dat import SliceMoms
from pp_h5dat import mkdirs_if_nexist
import pp_raw_ana
import h5py

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


def plot_save_proj_rms(slm, savepath):
    
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
    fig_sx.savefig( savepath + '/sigma_x_proj.eps',
                  format='eps')

    fig_sp = plt.figure()    
    plt.plot(t, np.sqrt(psq))
    ax = plt.gca()
    ax.set_xlabel(r'$\omega_p t$', fontsize=14) 
    ax.set_ylabel(r'$\sigma_{p_x}/m_e c$', fontsize=14)
    fig_sx.savefig( savepath + '/sigma_px_proj.eps',
                  format='eps')


    fig_xp = plt.figure()    
    plt.plot(t, xp)
    ax = plt.gca()
    ax.set_xlabel(r'$\omega_p t$', fontsize=14) 
    ax.set_ylabel(r'$ k_p x\,p_x/m_e c$', fontsize=14)    
    fig_xp.savefig( savepath + '/xpx_proj.eps',
                  format='eps')

    emittance = np.sqrt(xsq*psq-np.power(xp,2))
    fig_e = plt.figure()    
    plt.plot(t, emittance)
    ax = plt.gca()
    ax.set_xlabel(r'$\omega_p t$', fontsize=14) 
    ax.set_ylabel(r'$k_p \epsilon_x$', fontsize=14)    
    fig_e.savefig( savepath + '/emittance_proj.eps',
                  format='eps')

def plot_save_slice_centroids(slm, savepath):

    Xb0 = np.ones(slm.avgx2[0,:].shape)
    for i in range(0,len(slm.zeta_array)):
        if (slm.avgx2[0,i] != 0):
            Xb0[i] = slm.avgx2[0,i]


    Xb_norm = np.zeros( slm.avgx2.shape )
    zeta_hseed = 1.0
    idx_hseed = (np.abs(slm.zeta_array-zeta_hseed)).argmin()

    for i in range(0,len(slm.zeta_array)):
        if (slm.zeta_array[i] <= zeta_hseed):
            Xb_norm[:,i] = np.absolute( ( slm.avgx2[:,i] - slm.avgx2[:,idx_hseed])/Xb0[i] )

    fig1 = plt.figure()
    cax = plt.pcolormesh(   slm.zeta_array,
                            slm.time_array,
                            slm.avgx2,
                            cmap=cm.PuOr,
                            vmin=-np.amax(abs(slm.avgx2)), vmax=np.amax(abs(slm.avgx2)))
    cbar = fig1.colorbar(cax)
    cbar.ax.set_ylabel('$k_p X_b$')
    fig1.savefig( savepath + '/Xb_raw.png',
                  format='png',
                  dpi=600)

    fig2 = plt.figure()
    cax = plt.pcolormesh( Xb_norm )
    cbar = fig2.colorbar(cax)
    cbar.ax.set_ylabel('$|X_b/X_{b,0}|$')
    fig2.savefig( savepath + '/Xb.png',
                  format='png',
                  dpi=600)

    fig3 = plt.figure()
    plt.plot(slm.zeta_array, slm.avgx2[0,:])
    fig3.savefig( savepath + '/Xb0.png',
                  format='png')


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


def main():

    parser = ps_parseopts()

    args = parser.parse_args()

    file = args.path

    slm = SliceMoms(file)

    mkdirs_if_nexist(args.savepath)

    plot_save_slice_rms(slm, args.savepath)

    plot_save_slice_ene(slm, args.savepath)

    plot_save_proj_rms(slm, args.savepath)

    # plot_save_slice_centroids(slm, args.savepath)


if __name__ == "__main__":
    main()
