#!/usr/bin/env python3

import math
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.constants as constants
import plot_g3d


def add_parseargs(parser):
    parser.add_argument(  "--plasma-Z",
                          action='store',
                          dest="plasma_Z",
                          metavar="Z",
                          default=-1,
                          type=int,
                          help= """Atomic number of plasma species (Z=-1 for electrons).""")
    parser.add_argument(  "--plasma-A",
                          action='store',
                          dest="plasma_A",
                          metavar="A",
                          default=0,
                          type=int,
                          help= """Atomic number of plasma species (A=0 for electrons).""")
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
                          default=0,
                          type=float,
                          help= """Beam density""") 
    parser.add_argument(  "--beam-sigma_z",
                          action='store',
                          dest="beam_sigma_z",
                          metavar="sigma_z",
                          default=0,
                          type=float,
                          help= """Gaussian beam sigma_z""") 
    parser.add_argument(  "--beam-sigma_r",
                          action='store',
                          dest="beam_sigma_r",
                          metavar="sigma_r",
                          default=0,
                          type=float,
                          help= """Gaussian beam sigma_r""") 
    parser.add_argument(  "--beam-zeta_0",
                          action='store',
                          dest="beam_zeta_0",
                          metavar="zeta_0",
                          default=0.0,
                          type=float,
                          help= """Beam center location.""")                                                                     
    return parser                                                                

class Beam:
    def __init__(self, n=None, sigma_z = None, sigma_r = None, zeta_0 = 0.0):    
        if not any([n,sigma_z]) == None:
            self.n = n
            self.sigma_z = sigma_z
            self.sigma_r = sigma_r
            self.zeta_0 = zeta_0
        else:
            print('ERROR: Not all beam parameters specified!')

class Plasma:
    def __init__(self, n=1.0, A = 0, Z = -1):    
            self.n = n
            self.A = A
            self.Z = Z
            if self.A != 0:
                self.omega_p = np.sqrt(constants.m_e/(self.A * constants.m_p)) * abs(self.Z) * np.sqrt(self.n)
            else:
                self.omega_p = np.sqrt(self.n) * abs(self.Z)
            self.k_p = 1.0/self.omega_p

def lin_Ez_theo( plasma, beam, zeta_array ):

    Ez = -np.sign(plasma.Z) * beam.n/plasma.n * np.sqrt(2 * math.pi) * plasma.k_p * beam.sigma_z \
         * np.exp(-1 * (plasma.k_p * beam.sigma_z)**2/2) * np.cos( plasma.k_p * (zeta_array - beam.zeta_0) )

    return Ez


def cmp_plot_Ez(g3d_p, 
                plasma, 
                beam ):

    zeta_array = g3d_p.x_array

    Ez_theo = lin_Ez_theo( plasma, beam, zeta_array )

    fig = plt.figure()
    cax = plt.plot( g3d_p.x_array, g3d_p.line)
    cax = plt.plot( g3d_p.x_array, Ez_theo, '--')

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
    beam = Beam(n = args.beam_n, 
                sigma_z = args.beam_sigma_z, 
                sigma_r = args.beam_sigma_r,
                zeta_0 = args.beam_zeta_0 )
    return beam


def main():
    parser = plot_g3d.parser( ptype = 'line' )
    parser = add_parseargs( parser )

    args = parser.parse_args()
    
    beam = set_beam(args)
    plasma = set_plasma(args)


    flist = plot_g3d.gen_filelist( args )

    for file in flist:
        g3d_p = plot_g3d.G3d_plot_line(file, args)
        cmp_plot_Ez(g3d_p, plasma, beam)

if __name__ == "__main__":
    main()
