#!/usr/bin/env python3

import sys
import math
import argparse
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import scipy.constants as constants
import mpmath
import plot_g3d


def add_parseargs(plot_g3d_parser):
    parser = argparse.ArgumentParser(parents=[plot_g3d_parser],
                                     prog='plot_g3d_lin_vs_theo', 
                                     conflict_handler='resolve')
    parser.add_argument(  "--plasma-Z",
                          action='store',
                          dest="plasma_Z",
                          metavar="Z",
                          default=-1,
                          type=int,
                          help= """Atomic (charge) number of plasma species (Z=-1 for electrons).""")
    parser.add_argument(  "--plasma-A",
                          action='store',
                          dest="plasma_A",
                          metavar="A",
                          default=0,
                          type=int,
                          help= """Mass number of plasma species (A=0 for electrons).""")
    parser.add_argument(  "--plasma-n",
                          action='store',
                          dest="plasma_n",
                          metavar="np",
                          default=1.0,
                          type=float,
                          help= """Plasma density""") 
    parser.add_argument(  "--beam-n",
                          action='store',
                          dest="beam_n",
                          metavar="nb",
                          default=0.0,
                          type=float,
                          required=True,
                          help= """Beam density""") 
    parser.add_argument(  "--beam-sigma_z",
                          action='store',
                          dest="beam_sigma_z",
                          metavar="sigma_z",
                          default=0.0,
                          type=float,
                          required=True,
                          help= """Gaussian beam sigma_z""") 
    parser.add_argument(  "--beam-sigma_r",
                          action='store',
                          dest="beam_sigma_r",
                          metavar="sigma_r",
                          default=0.0,
                          type=float,
                          help= """Gaussian beam sigma_r""") 
    parser.add_argument(  "--beam-sigma_xy",
                          action='store',
                          dest="beam_sigma_xy",
                          metavar="sigma_xy",
                          default=0.0,
                          type=float,
                          help= """Gaussian beam sigma_x = sigma_y""")     
    parser.add_argument(  "--beam-zeta_0",
                          action='store',
                          dest="beam_zeta_0",
                          metavar="zeta_0",
                          default=0.0,
                          type=float,
                          help= """Beam center location.""")
    parser.add_argument(  "--beam-Z",
                          action='store',
                          dest="beam_Z",
                          metavar="Z",
                          default=-1,
                          type=int,
                          help= """Beam atomic (charge) number.""")
    parser.add_argument(  "-s", "--save-path",
                          action="store",
                          dest="savepath",
                          metavar="PATH",
                          default=plot_g3d.parsedefs.save.path + '/g3d-line_vs_theo',
                          help = """Path to which generated files will be saved.
                              (Default: './')""")
    parser.add_argument(  "-f", "--format",
                          action='store',
                          dest="file_format",
                          metavar="FORMAT",
                          choices=[ 'png',
                                    'pdf',
                                    'eps',],
                          default='eps',
                          help= """Format of output file (Default: eps).""")

    parser.set_defaults(lout_idx=None)
    parser.set_defaults(loutax='z')

    return parser                                                                

class Beam:
    def __init__(self, n=0.0, sigma_z = 0.0, sigma_r = 0.0, zeta_0 = 0.0, Z=-1):    
        if not (n == 0.0 or sigma_z == 0.0):
            self.n = n
            self.sigma_z = sigma_z
            self.sigma_r = sigma_r
            self.zeta_0 = zeta_0
            self.Z = Z
        else:
            print('ERROR: Beam density and beam sigma_z must be specified!')
            sys.exit()

class Plasma:
    def __init__(self, 
                 n = 1.0, 
                 A = 0, 
                 Z = -1):    
            self.n = n
            self.A = A
            self.Z = Z

            if self.A != 0:
                print('Plasma A = %d' % A)
                self.omega_p = np.sqrt(constants.m_e/(self.A * constants.m_p)) * abs(self.Z) * np.sqrt(self.n)
            else:
                self.omega_p = np.sqrt(self.n) * abs(self.Z)
            self.k_p = self.omega_p
            self.lambda_p = 2 * math.pi / self.k_p
            print('k_p = %f' % self.k_p)

def lin_Ez_theo_long_wave( plasma, beam, zeta_array ):

    Ez = -np.sign(beam.Z) * np.sqrt(2 * math.pi) * beam.n/plasma.n * plasma.k_p * beam.sigma_z \
         * np.exp(-1 * (plasma.k_p * beam.sigma_z)**2/2) * np.cos( plasma.k_p * (zeta_array - beam.zeta_0) )

    return Ez


def lin_Ez_theo_sigma_r( plasma, beam, zeta_array ): 

    Ez = -np.sign(beam.Z) * np.sqrt(2 * math.pi) * beam.n/plasma.n * plasma.k_p*beam.sigma_z \
         * np.exp(-1 * (plasma.k_p*beam.sigma_z)**2/2) * np.cos( plasma.k_p * (zeta_array - beam.zeta_0) ) \
         * (plasma.k_p*beam.sigma_r)**2/2 * np.exp((plasma.k_p*beam.sigma_r)**2/2) * mpmath.gammainc(0.0,(plasma.k_p*beam.sigma_r)**2/2)

    return Ez


def cmp_plot_Ez(g3d_p, 
                plasma, 
                beam ):

    zeta_array = g3d_p.x_array

    if beam.sigma_r == 0.0 or beam.sigma_r == None: 
        Ez_theo = lin_Ez_theo_long_wave( plasma, beam, zeta_array )
        print('lin_Ez_theo_long_wave')
    else:
        Ez_theo = lin_Ez_theo_sigma_r( plasma, beam, zeta_array )

    print('Ratio: %f' % (np.max(Ez_theo)/np.max(g3d_p.line) ) )

    fig = plt.figure()
    ax = plt.plot( g3d_p.x_array, g3d_p.line )
    ax = plt.plot( g3d_p.x_array, Ez_theo, '--' )
    ax = plt.gca()
    ax.set_ylabel(g3d_p.ylabel, fontsize=14)
    ax.set_xlabel(g3d_p.xlabel, fontsize=14)

    if not (-3.0 < math.log(np.max(abs(g3d_p.line)),10) < 3.0):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15)  

    g3d_p.mkdirs_if_nexist()

    saveformat = g3d_p.args.file_format
    filesuffix = '_%06.f' % (np.floor(g3d_p.g3d.time))

    fileprefix = g3d_p.g3d.name

    savename = fileprefix + filesuffix + '.' + saveformat

    if saveformat=='png':
        fig.savefig(  g3d_p.args.savepath + '/' + savename,
                  format=saveformat,
                  dpi=600)
    else:
        fig.savefig(  g3d_p.args.savepath + '/' + savename,
                      format=saveformat)
    if g3d_p.args.verbose: print('Saved "' + savename + '" at: ' + g3d_p.args.savepath)

    if g3d_p.args.ifshow: plt.show()
    plt.close(fig)


def set_plasma( args ):
    plasma = Plasma(n = args.plasma_n, 
                    A = args.plasma_A, 
                    Z = args.plasma_Z )
    return plasma


def set_beam( args ):
    if args.beam_sigma_r != 0.0 and args.beam_sigma_xy == 0.0:
        beam = Beam(n = args.beam_n, 
                    sigma_z = args.beam_sigma_z, 
                    sigma_r = args.beam_sigma_r,
                    zeta_0 = args.beam_zeta_0,
                    Z = args.beam_Z )
    elif args.beam_sigma_r == 0.0 and args.beam_sigma_xy != 0.0:
        beam = Beam(n = args.beam_n, 
                    sigma_z = args.beam_sigma_z, 
                    sigma_r = args.beam_sigma_xy * np.sqrt(2), # DOUBLE CHECK THIS!!!
                    zeta_0 = args.beam_zeta_0,
                    Z = args.beam_Z )
    else:
        print('ERROR: "beam_sigma_r" can''t be set in conjunction with "beam_sigma_xy"!')
        sys.exit(1)     
    return beam


def main():
    plot_g3d_parser = plot_g3d.g3d_parser()
    parser = add_parseargs( plot_g3d_parser )
    args = parser.parse_args()
    
    beam = set_beam(args)
    plasma = set_plasma(args)


    flist = plot_g3d.gen_filelist( args )

    for file in flist:
        g3d_p = plot_g3d.G3d_plot_line(file, args)
        cmp_plot_Ez(g3d_p, plasma, beam)

if __name__ == "__main__":
    main()
