#!/usr/bin/env python3

import sys
import math
import argparse
import numpy as np
import itertools
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter
import scipy.constants as constants
import mpmath

mypath = os.path.dirname(os.path.abspath( __file__ ))
incpath = os.path.split(mypath)[0] + '/inc'
sys.path.append(incpath)
from pp_h5dat import Grid3d
from pp_h5dat import H5FList
from pp_h5dat import H5Plot
from pp_h5dat import mkdirs_if_nexist
from pp_plt_tools import saveas_eps_pdf

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
        prefix = 'g3d_name'
        path = './plots'


def flip(items, ncol):
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])

def g3d_lvst_parser():

    desc = """This is the picpy line vs. theory plotting tool."""
    # Line vs-theo plot arguments
    parser = argparse.ArgumentParser( add_help=False,
                                      description=desc)
    parser.add_argument(  'path',
                          metavar = 'PATH',
                          nargs = '*',
                          help = 'Path to grid file.')
    parser.add_argument(  "-v", "--verbose",
                          dest = "verbose",
                          action="store_true",
                          default=True,
                          help = "Print info (Default).")
    parser.add_argument(  "-q", "--quiet",
                          dest = "verbose",
                          action = "store_false",
                          help = "Don't print info.")    
    parser.add_argument(  "-s", "--save-path",
                          action="store",
                          dest="savepath",
                          metavar="PATH",
                          default=parsedefs.save.path + '/g3d-line-vs-theo',
                          help = """Path to which generated files will be saved.
                              (Default: %(default)s)""")
    parser.add_argument(  "-f", "--format",
                          action='store',
                          dest="file_format",
                          metavar="FORMAT",
                          choices=[ 'png',
                                    'pdf',
                                    'eps',],
                          default='eps',
                          help= """Format of output file (Default: %(default)s).""")
    parser.add_argument(  "--h5",
                          dest = "h5plot",
                          action="store_true",
                          default=True,
                          help = "Save plot as hdf5 file (Default: %(default)s).")    
    parser.add_argument(  "--ion-Z",
                          action='store',
                          dest="ion_Z",
                          metavar="Z",
                          default=-1,
                          type=int,
                          help= """Atomic (charge) number of ion species (Z=-1 for electrons).""")
    parser.add_argument(  "--ion-A",
                          action='store',
                          dest="ion_A",
                          metavar="A",
                          default=0,
                          type=int,
                          help= """Mass number of ion species (A=0 for electrons).""")
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
    return parser


def g3d_lvst_Ez_subparser(subparsers, g3d_lvst_parent):
    parser = subparsers.add_parser( "Ez", parents=[g3d_lvst_parent],
                                    conflict_handler='resolve',
                                    help="Ez longitudinal lineout")
    parser.add_argument( "-l", "--lineout-axis",
                      action='store',
                      dest="loutax",
                      metavar="LOUTAX",
                      choices=[ 'x', 'y', 'z'],
                      default='z',
                      help= """Axis along which lineout is generated (Default: z).""")
    parser.add_argument(  '--lineout-indices',
                          help='Indices for which lineout is taken.',
                          action='store',
                          dest="lout_idx",
                          metavar=('IDX1', 'IDX2'),
                          type=int,
                          nargs=2,
                          default=None)
    parser.add_argument(  "--beam-sigma_z",
                          action='store',
                          dest="beam_sigma_z",
                          metavar="sigma_z",
                          default=0.0,
                          type=float,
                          help= """Gaussian beam sigma_z""")                                                      
    return parser


def g3d_lvst_Wr_subparser(subparsers, g3d_lvst_parent):
    parser = subparsers.add_parser( "Wr", parents=[g3d_lvst_parent],
                                    conflict_handler='resolve',
                                    help="Wr transverse lineout")
    parser.add_argument( "-l", "--lineout-axis",
                      action='store',
                      dest="loutax",
                      metavar="LOUTAX",
                      choices=[ 'x', 'y', 'z'],
                      default='x',
                      help= """Axis along which lineout is generated (Default: x).""")
    parser.add_argument(  "--zeta-pos",
                      action='store',
                      dest="zeta_pos",
                      metavar="LOUT-ZETA-POS",
                      type=float,                      
                      default=None,
                      nargs='+',
                      help= """Zeta-position at which lineout is generated (Default: 0.0).""")
    parser.add_argument(  '--lineout-indices',
                          help='Indices for which lineout is taken.',
                          action='store',
                          dest="lout_idx",
                          metavar=('IDX0', 'IDX1'),
                          type=int,
                          nargs=2,
                          default=None)
    parser.add_argument(  "--rmax",
                      action='store',
                      dest="rmax",
                      metavar="RMAX",
                      type=float,                      
                      default=None,
                      help= """Maximum radius for plotting.""")
    parser.add_argument(  "--rel-to-hom",
                          dest = "rel_to_hom",
                          action="store_true",
                          default=False,
                          help = "Plot Wr-r/2 (Default: %(default)s).")                        
    return parser

class Beam:
    def __init__(self, n=0.0, sigma_z = 0.0, sigma_xy = 0.0, zeta_0 = 0.0, Z=-1, L=None, n_head=None, n_tail=None):    
        if not (n == 0.0):
            self.n = n
            self.sigma_z = sigma_z
            self.sigma_xy = sigma_xy
            self.zeta_0 = zeta_0
            self.Z = Z
            self.L = L

        else:
            print('ERROR: Beam density must be specified!')
            sys.exit()

class Plasma:
    def __init__(self, 
                 n = 1.0, 
                 A = 0, 
                 Z = -1):    
            self.n = n
            self.A = A
            self.Z = Z
            # if self.A != 0:
            #     print('Ion A = %d' % A)
            #     print('Ion Z = %d' % Z)
            #     self.omega_p = np.sqrt(constants.m_e/(self.A * constants.m_p)) * abs(self.Z) * np.sqrt(self.n)
            # else:
            #     self.omega_p = np.sqrt(self.n) * abs(self.Z)
            # self.k_p = self.omega_p
            # self.lambda_p = 2 * math.pi / self.k_p
            # print('k_p = %f' % self.k_p)

def lin_Ez_theo_long_pwave( plasma, beam, zeta_array ):

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
        Ez_theo = lin_Ez_theo_long_pwave( plasma, beam, zeta_array )
        print('using: lin_Ez_theo_long_pwave')
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
    filesuffix = '_%06.f' % (np.floor(g3d_p.g3d.get_time()))

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



def lin_Wr_theo( plasma, beam, x_array, zeta_pos ): 

    def H(x):
        return (1-np.exp(-x))/x

    if plasma.A == 0:
        Wr = x_array/2
    else:
        Wr = x_array/2 * ( 1.0 + plasma.Z 
                               * constants.m_e/(plasma.A * constants.m_p)
                               * beam.n 
                               * zeta_pos**2/2 
                               * H( x_array**2/( 2 * beam.sigma_xy**2 ) ) )

    return Wr

def cmp_plot_Wr(args,
                g3d, 
                plasma, 
                beam ):

    if args.zeta_pos == None:
        zeta_pos_list = [0.0]
    else:
        zeta_pos_list = args.zeta_pos

    x_array = g3d.get_x_arr(1)
       
    Nsimlines = len(zeta_pos_list)
    Nzeta = g3d.get_nx(1)

    Wr_sim = np.zeros((Nzeta, Nsimlines))
    Wr_theo = np.zeros((Nzeta, Nsimlines))

    for i in range(0, Nsimlines):
        zeta_pos = zeta_pos_list[i]
        Wr_sim[:,i] = g3d.read(x0 = zeta_pos, x2 = 0.0)
        Wr_theo[:,i] = lin_Wr_theo( plasma = plasma, 
                               beam = beam,
                               x_array = x_array,
                               zeta_pos = zeta_pos)

    cmap = cm.tab10


    fig = plt.figure(num=None, 
                     figsize=(9, 7), 
                     dpi=80, 
                     facecolor='w', 
                     edgecolor='k') 
    for i in range(0, Nsimlines):
        zeta_pos = zeta_pos_list[i]
        if args.rel_to_hom:
            Wr_sim_vals = Wr_sim[:,i] - x_array*0.5
            Wr_theo_vals = Wr_theo[:,i] - x_array*0.5
            ylab = r'$W_r/E_0 - r/2$'
        else:
            Wr_sim_vals = Wr_sim[:,i]
            Wr_theo_vals = Wr_theo[:,i]
            ylab = r'$W_r/E_0$'            

        label_sim = r'PIC: $k_p \zeta = %0.1f$' %  zeta_pos
        ax_sim = plt.plot( x_array, Wr_sim_vals,
                           linestyle='-',
                           color=cmap(i),
                           label=label_sim)    
        
        label_theo = r'Theo: $k_p \zeta = %0.1f$' %  zeta_pos
        ax_theo = plt.plot( x_array, Wr_theo_vals, 
                            linestyle ='--',
                            color=cmap(i),
                            label=label_theo)
    if not args.rel_to_hom:
        label_half = r'$k_p r/2$'
        ax_half = plt.plot( x_array, x_array*0.5, 
                            linestyle ='-',
                            label=label_half,
                            color=[0.5, 0.5, 0.5])

    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(flip(handles, 2), flip(labels, 2), ncol=2)

    if args.rmax != None:
        xmax = args.rmax
        ymax = 1.2 * np.amax(Wr_sim_vals[np.logical_and(x_array<xmax, x_array>0.0)])
    else:
        max_idx = np.where(Wr_sim_vals==np.amax(Wr_sim_vals))[0][0]
        xmax = x_array[max_idx]
        ymax = 1.2 * np.amax(Wr_sim_vals)
    ax.set_xlim([0,xmax])
    ax.set_ylim([0,ymax])            
    ax.set_ylabel(ylab, fontsize=14)
    ax.set_xlabel(r'$k_p r$', fontsize=14)

    if not (-3.0 < math.log(np.max(abs(Wr_sim)),10) < 3.0):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15)  

    filesuffix = '_t_%06.f' % (np.floor(g3d.get_time()))
    fileprefix = 'Wr_' + 'zeta_%0.1f' % zeta_pos
    saveformat = 'pdf'
    savepath = args.savepath

    savename = fileprefix + filesuffix + '.' + saveformat

    saveas_eps_pdf(fig, savepath, savename, h5plot=args.h5plot, fformat=saveformat)
    
    plt.close(fig)


def set_plasma( args ):
    plasma = Plasma(n = args.plasma_n, 
                    A = args.ion_A, 
                    Z = args.ion_Z )
    return plasma


def set_beam( args ):

    if args.subparser_name == 'Ez':
        sigma_z = args.beam_sigma_z
    else:
        sigma_z = 0

    if args.beam_sigma_r != 0.0 and args.beam_sigma_xy == 0.0:
        sigma_xy = args.beam_sigma_r/np.sqrt(2.0)

    elif args.beam_sigma_r == 0.0 and args.beam_sigma_xy != 0.0:
        sigma_xy = args.beam_sigma_xy
        
    elif args.beam_sigma_r == 0.0 and args.beam_sigma_xy == 0.0:
        print('ERROR: Either "beam_sigma_r" or "beam_sigma_xy" must be set!')
        sys.exit(1)  
    else:
        print('ERROR: "beam_sigma_r" can''t be used in conjunction with "beam_sigma_xy"!')
        sys.exit(1)
    beam = Beam(n = args.beam_n, 
                sigma_z = sigma_z, 
                sigma_xy = sigma_xy,
                zeta_0 = args.beam_zeta_0,
                Z = args.beam_Z )          
    return beam



def lvst_Ez(args):  
    beam = set_beam(args)
    plasma = set_plasma(args)

    flist = plot_g3d.gen_filelist( args )

    for file in flist:
        g3d_p = plot_g3d.G3d_plot_line(file, args)
        if (g3d_p.g3d.name == 'Ez'):  
            cmp_plot_Ez(g3d_p, plasma, beam)
        else:
            print('ERROR: Dataset is no Ez dataset!')
            sys.exit(1)                 


def lvst_Wr(args):  
    beam = set_beam(args)
    plasma = set_plasma(args)

    h5flist = H5FList(args.path, 'g3d')
    flist = h5flist.get()

    for file in flist:
        g3d = Grid3d(file)
        if (g3d.name == 'ExmBy' or g3d.name == 'EypBx'):  
            cmp_plot_Wr(args, g3d, plasma, beam)
        else:
            print('ERROR: Dataset is no ExmBy or EypBx dataset!')
            sys.exit(1)    

def main():

    parser = argparse.ArgumentParser()
    g3d_subparsers = parser.add_subparsers(title="plot-type",
                                           dest="subparser_name")

    g3d_lvst_pp = g3d_lvst_parser()

    g3d_lvst_Ez_sp = g3d_lvst_Ez_subparser( subparsers=g3d_subparsers,
                                            g3d_lvst_parent=g3d_lvst_pp)
    g3d_lvst_Ez_sp.set_defaults(func=lvst_Ez)

    g3d_lvst_Wr_sp = g3d_lvst_Wr_subparser( subparsers=g3d_subparsers,
                                            g3d_lvst_parent=g3d_lvst_pp)
    g3d_lvst_Wr_sp.set_defaults(func=lvst_Wr)

    args = parser.parse_args()
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    args.func(args)


if __name__ == "__main__":
    main()
