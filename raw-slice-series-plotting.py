#!/usr/bin/env python3
# This script may be executed like this:
# nohup ./raw-slice-series-plotting.py <DATA>/ 1> rss.out 2> rss.err &

import numpy as np
import os
import argparse
from optparse import OptionParser
#from optparse import OptionGroup
import math
import sys
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import picdefs
from h5dat import SliceMoms
import ps_ana
import h5py

# Parse defaults/definitions
class parsedefs:
    class file_format:
        png = 'png'
        eps = 'eps'
        pdf = 'pdf'
    class save_prefix:
        name = 'g3d_name'


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
                        help = "Print info (Default).")
    parser.add_argument(  "-q", "--quiet",
                        action="store_false",
                        dest="verbose",
                        help = "Don't print info.")
    parser.add_argument(  "-s", "--save-path",
                        dest="savepath",
                        metavar="PATH",
                        default='./',
                        help = """Path to which generated files will be saved.
                              (Default: './')""")
    parser.add_argument(  "-n", "--name-prefix",
                        dest="save_prefix",
                        metavar="NAME",
                        default=parsedefs.save_prefix.name,
                        help = """Define customized prefix of output filename.""")
    parser.add_argument(  "-c", "--code",
                        action='store',
                        dest="piccode",
                        metavar="CODE",
                        choices = [picdefs.code.hipace, picdefs.code.osiris,],
                        default = picdefs.code.hipace,
                        help= "PIC code which was used to generate files (Default: " +
                              picdefs.code.hipace + ").")
    parser.add_argument(  "-d", "--dim",
                        action='store',
                        dest="dimensionality",
                        metavar="DIM",
                        choices=[1, 2, 3,],
                        default=3,
                        help= """Dimensionality of PIC simulation
                              (Default: 3).""")
    parser.add_argument(  "-N", "--number-of-files",
                        action='store',
                        dest="Nfiles",
                        metavar="NFILES",
                        default=0,
                        help= """Number of files to analyze.""")

    return parser



def main():

    parser = ps_parseopts()

    args = parser.parse_args()

    file = args.path

    slm = SliceMoms(file)

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
    cax = plt.pcolormesh( slm.avgx2 )
    cbar = fig1.colorbar(cax)
    cbar.ax.set_ylabel('$X_b$')
    fig1.savefig(  'Xb_raw.png',
                  format='png')

    fig2 = plt.figure()
    cax = plt.pcolormesh( Xb_norm )
    cbar = fig2.colorbar(cax)
    cbar.ax.set_ylabel('$|X_b/X_{b,0}|$')
    fig2.savefig(  'Xb.png',
                  format='png')

    fig3 = plt.figure()
    plt.plot(slm.zeta_array, slm.avgx2[0,:])
    fig3.savefig(  './Xb0.png',
                  format='png')


    fig4 = plt.figure()
    cax = plt.pcolormesh( slm.zeta_array,
                          slm.time_array,
                          np.sqrt( np.absolute( slm.avgx2sq ) ) )
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$\omega_p t$', fontsize=14)    
    cbar = fig4.colorbar( cax )
    cbar.ax.set_ylabel(r'$k_p \sigma_x$', fontsize=14)
    fig4.savefig(  'sigma_x.png',
                  format='png')

    emittance = np.sqrt( np.multiply(slm.avgx2sq, slm.avgp2sq) 
                         - np.power(slm.avgx2p2,2))

    fig5 = plt.figure()
    cax = plt.pcolormesh( slm.zeta_array,
                          slm.time_array,
                          emittance)
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$\omega_p t$', fontsize=14)    
    cbar = fig5.colorbar( cax )
    cbar.ax.set_ylabel(r'$k_p \epsilon_x$', fontsize=14)
    fig5.savefig(  'emittance_x.png',
                  format='png')


    # fig5 = plt.figure()
    # cax = plt.plot( slm.time_array, np.sqrt( np.absolute( slm.avgx2sq[:,80] )) )
    # cbar.ax.set_ylabel('$\sigma_x$')
    # fig5.savefig(  'sigma_x_slice.png',
    #                 format='png')

if __name__ == "__main__":
    main()
