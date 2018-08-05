#!/usr/bin/env python3
import csv
import os
import glob
import numpy as np
import argparse
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')
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

from decimal import *

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

    parser.add_argument(  "--box",
                          dest = "rmax",
                         # action="store_true",
                          type = float,
                          default=0,
                          help = "max transverse box size (Default: %(default)s).")
                          
    parser.add_argument(    "--nb",
                            dest = "nb",
                        #    action="store_true",
                            type = float,
                            default=0,
                            help = "bunch density (Default: %(default)s).")
                    
    parser.add_argument(  "--nx",
                          dest = "nx",
                         # action="store_true",
                          type = float,
                          default=0,
                          help = "number of transverse gridpoints in hipace input script (Default: %(default)s).")
                          
    parser.add_argument(  "--ni",
                        dest = "ni",
                        #action="store_true",
                        type = float,
                        default=1,
                        help = "ion density (Default: %(default)s).")
                        
    parser.add_argument(  "--rbunch",
                          dest = "rbunch",
                         # action="store_true",
                          type = float,
                          default=0,
                          help = "Radius or Sigma_r of the bunch as in hipace input script (Default: %(default)s).")

    parser.add_argument(  "--lbunch",
                        dest = "lbunch",
                        #action="store_true",
                        type = float,
                        default=2,
                        help = "bunch length (Default: %(default)s).")
    parser.add_argument(  "--zeta-pos",
                        dest = "zeta_pos",
                        #action="store_true",
                        type = float,
                        default= -0.2,
                        help = "bunch length (Default: %(default)s).")
    parser.add_argument(  "--I_beam",
                        dest = "I_beam",
                        #action="store_true",
                        type = float,
                        default= 0.5,
                        help = "bunch current in I_b/I_A (Default: %(default)s).")

    parser.add_argument(  "--ionization",
                        dest = "ionization",
                        #action="store_true",
                         action="store_true",
                         default=False,
                         help = "Save and plot ionization probability as hdf5 file (Default: %(default)s).")
    parser.add_argument(  "--canbeam",
                        dest = "canbeam",
                        #action="store_true",
                         action="store_true",
                         default=False,
                         help = "Save and plot ionization probability as hdf5 file (Default: %(default)s).")


    return parser


def display_cmap(cmap):
    print("display_cmap disabled")
    # plt.imshow(np.linspace(0, 100, 256)[None, :],  aspect=25,    interpolation='nearest', cmap=cmap) 
    # plt.axis('off')
    # plt.show()

def calc_ion_rate( elec_field,
                   Z,
                   ion_pot,
                   ion_pot_zero,
                   l,
                   abs_m): 
    #print('Z = %f' % Z)
    #print('ion_pot = %f' % ion_pot)
    #print('U_h = %f' % HYDROGEN_IONIZATION_ENERGIE_EV)
    omega_alpha               =    4.13e16; # [s^-1]
    E_alpha                   =    5.1e11; 
    HYDROGEN_IONIZATION_ENERGIE_EV = 13.659843449
    n_eff = Z * np.sqrt(HYDROGEN_IONIZATION_ENERGIE_EV/ion_pot)
    l_eff = np.sqrt(HYDROGEN_IONIZATION_ENERGIE_EV/ion_pot_zero) - 1
    C_2 =  2**(2*n_eff)/(n_eff *scipy.special.gamma(n_eff + l_eff + 1) * scipy.special.gamma(n_eff - l_eff) ) ## approx 1/(2*np.pi*n_eff) * (2/n_eff)**(2*n_eff)#
    quant_prefactor = ((2*l + 1)*scipy.math.factorial(l+abs_m)) / ( 2**abs_m * scipy.math.factorial(abs_m) * scipy.math.factorial(l - abs_m))
                              
    value = 0
    # print('n_eff = %f' % n_eff)
    # print('C_2 = %f' % C_2)
    # print('Quant_factor = %f' % quant_prefactor)
    value1 = 0;
    value2 = 0;
    value3 = 0;
    elec_field = np.abs(elec_field)
    if (elec_field/E_alpha > 1e-3):
        # print('hier bin ich richtig')
        value1 = omega_alpha * quant_prefactor * C_2 * ion_pot/ (2* HYDROGEN_IONIZATION_ENERGIE_EV)
        value2 = (2*E_alpha/elec_field * (ion_pot/ HYDROGEN_IONIZATION_ENERGIE_EV)**1.5)**(2*n_eff -1) 
        value3 = np.exp(-2.0/3.0 * E_alpha/elec_field * (ion_pot/ HYDROGEN_IONIZATION_ENERGIE_EV)**1.5 )
        # print('value 1 : %f' %value1)
        # print('value 2 : %f' %value2)
        # print('value 3 : %f' %value3)
    return value1*value2*value3 


def plot_hipace_Ex(zeta_pos):
    ExmBy_path = './DATA/field_ExmBy_000000.0.h5'
    By_path = './DATA/field_By_000000.0.h5'
    
    ExmBy_g3d2 = Grid3d(ExmBy_path)
    By_g3d2 = Grid3d(By_path)
    
    ExmBy = np.transpose(ExmBy_g3d2.read(x2=0.0))
    By = np.transpose(By_g3d2.read(x2=0.0))  

    Ex = ExmBy + By
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    idx =np.abs(By_g3d2.get_zeta_arr() - zeta_pos).argmin()
    
    return By_g3d2.get_x_arr(2), Ex[:,idx]
    
def gauss_E_field(r, sigma, I_b):
    return -2*I_b*1/r * (1-np.exp(-r**2/(2*sigma**2)))

def calc_transversal_probability_density(r_max, nx, nb, ni, rbunch, lbunch, zeta_pos, I_beam):
    ELECTRON_CHARGE_IN_COUL   =    1.60217657e-19
    ELECTRON_MASS_IN_KG       =    9.10938291e-31
    VAC_PERMIT_FARAD_PER_M    =    8.854187817e-12
    SPEED_OF_LIGHT_IN_M_PER_S =    299792458
    omega_alpha               =    4.13e16; # [s^-1]
    E_alpha                   =    5.1e11; # [V/m]
    elec_density              =    1e+23 # [1/m^3]
    omega_p = np.sqrt(4* np.pi * elec_density * (ELECTRON_CHARGE_IN_COUL**2)/ (VAC_PERMIT_FARAD_PER_M * ELECTRON_MASS_IN_KG));
    E_0 = omega_p * ELECTRON_MASS_IN_KG * SPEED_OF_LIGHT_IN_M_PER_S / ELECTRON_CHARGE_IN_COUL;     
    #Use zeta array in order to calculate delta t
    
    savepath = './plots/g3d-line/'
    mkdirs_if_nexist(savepath)
    
    bunch_array = np.arange(0, lbunch+lbunch/10, lbunch/10)
    bunch_array_sum = np.cumsum(bunch_array)

    deltat = np.abs(bunch_array_sum) / omega_p

    vcalc_ion_rate = np.vectorize(calc_ion_rate)

    r_array = np.arange(0, r_max+r_max/nx, r_max/nx)

    
    c1prime = -1/2* nb 
    c3prime = -1/2*rbunch**2*(nb)

    
    savename = 'Ex_analytic'
    
    E = np.zeros(shape = np.shape(r_array))
    for i in range(len(r_array)):
        if args.canbeam:
            if(r_array[i] <= rbunch):
                E[i] = np.real(c1prime) * r_array[i] 
            else:
                E[i] = np.real(c3prime) / r_array[i]
        else:
            E[i] = gauss_E_field(r_array[i], rbunch, I_beam)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.plot(np.append(-r_array[::-1], r_array), np.append(-(E[::-1]) ,E[:]))
    x, Ex =plot_hipace_Ex(zeta_pos) # input = zeta pos 
    ax.plot(x, Ex)
        # x, Ex =plot_hipace_Ex(-0.2) # input = zeta pos 
        # ax.plot(x, Ex, label=r'$\zeta = -0.2$', 'r')
        # x, Ex =plot_hipace_Ex(-1.0) # input = zeta pos 
        # ax.plot(x, Ex, label=r'$\zeta = -1.0$', 'b')
        # x, Ex =plot_hipace_Ex(-2.0) # input = zeta pos 
        # ax.plot(x, Ex, label=r'$\zeta = -2.0$', 'g')
    
    
    
    h5lp = H5Plot()
    h5lp.inherit_matplotlib_line_plots(ax)
    h5lp.write(savepath + '/' + savename + '.h5')
    #plt.show()
    #use zeta array for shaping the elctric field to get the same spacing etc.
    #print(Earray)

####### CALCULATING IONIZATION   
    if args.ionization:
        ionization_rate =np.zeros(shape=np.shape(r_array))
        ionization_rate[:] = vcalc_ion_rate(E*E_0, 1, 13.659843449,13.659843449, 0,0 ) # complete formulae for hydrogen
        #print(ionization_rate[1])
        
        
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        #print('Deltat is %0.3e' %deltat)
        #computing the ionization probability
        
        
        savename = 'ion_probability_test'
        ion_probability = np.zeros(shape=(len(bunch_array) -1, len(ionization_rate)))
        for i in range(1,2): #len(bunch_array)):
            ion_probability[i-1,:] = 1.0 - np.exp(-ionization_rate[:] * deltat[i] )
            ax2.plot(np.append(-r_array[::-1], r_array), np.append((ion_probability[i-1,::-1]) ,ion_probability[i-1,:]), label=('prob at lb: ' + str(np.around(bunch_array[i], decimals = 1) ) ) )
            #ax.legend()
        
        
        h5lp = H5Plot()
        h5lp.inherit_matplotlib_line_plots(ax2)
        h5lp.write(savepath + '/' + savename + '.h5')
        # #plt.show()





def main():
    
    
    
    parser = binSlab_parser() 
    args = parser.parse_args()

    
    
    

    calc_transversal_probability_density(args.rmax, args.nx/2, args.nb, args.ni, args.rbunch, args.lbunch, args.zeta_pos, args.I_beam) #parameters from hipace.c ##(r_max, nx, nb, ni, rbunch, lbunch)


if __name__ == "__main__":
    main()  