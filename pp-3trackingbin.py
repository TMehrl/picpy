#!/usr/bin/env python3
import csv
import os
import glob
import numpy as np
import argparse
import matplotlib
matplotlib.use('TkAgg')
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
    parser.add_argument(  "--2Dproj",
                          dest = "twodproj",
                          default = False,
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


def erf_approx(x):
    a = 0.147
    
    value = np.sign(x)*np.sqrt( 1-np.exp( -x**2 * (4/np.pi + a*x**2)/(1+ a*x**2 )  )  )
    return value

def ierf_approx(x):
    a = 0.147
    value = np.sign(x)*np.sqrt( np.sqrt( (2/(np.pi*a) + np.log(1-x**2)/2)**2 - np.log(1-x**2)/a  ) - ( 2/(np.pi*a) + np.log(1-x**2)/2 )  )
    return value

def calc_ierfi(x):
    
    def f(y):
        return x - scipy.special.erfi(y)
    y = scipy.optimize.bisect(f, -100000,100000)
    return y

def calc_r(zeta, c, int_const_1, int_const_2):
    r = 1j * np.exp( -1j * np.sqrt(c)* (int_const_2 + zeta) ) * \
        ( np.exp( 2j*np.sqrt(c)*int_const_2 ) + 8*np.exp( 2j*np.sqrt(c)*zeta )* \
        ( 3*int_const_1 - np.log(2))  ) / ( 4*np.sqrt(3*c) )
    
    return r
    
def calc_r2_vp(zeta, c, int_const_1, int_const_2): #here I assumed v = - blabla because then you end up with this form
    r = 1j * np.exp( -1j * np.sqrt(c)* (int_const_2 + zeta) ) * \
        ( np.exp( 2j*np.sqrt(c)*int_const_2 ) - 8*np.exp( 2j*np.sqrt(c)*zeta )* \
        ( int_const_1 - 1))   / ( 4*np.sqrt(c) )
    
    return r

def calc_r2_vm(zeta, c, int_const_1, int_const_2): #here I assumed v = - blabla because then you end up with this form
    r = 1j * np.exp( -1j * np.sqrt(c)* (int_const_2 + zeta) ) * \
        ( np.exp( 2j*np.sqrt(c)*zeta ) - 8*np.exp( 2j*np.sqrt(c)*int_const_2 )* \
        ( int_const_1 - 1))   / ( 4*np.sqrt(c) )
    
    return r
    
    
def calc_r3_vm(zeta, c, int_const_1, int_const_2): #here I assumed v = - blabla because then you end up with this form
    
    if np.real(c)<0:
        # IERFI_APPROX = calc_ierfi( np.real( np.sqrt(np.abs(c))* np.exp( (1-int_const_1)/c ) * np.sqrt(2/np.pi) * (-int_const_2 + zeta) ) )
        # #print('IERFI_APPROX:  %f ' %IERFI_APPROX)
        # r = np.exp( (-1 + int_const_1 +c* IERFI_APPROX**2)/c )
        r = np.exp( (-1 + int_const_1 -c* ierf_approx(( np.sqrt(c)* np.exp( (1-int_const_1)/c ) * np.sqrt(2/np.pi) * (-int_const_2 + zeta) ))**2)/c )
    else:
        #print('INPUT IN ERFINV : %f imag %f' %(np.real(( np.sqrt(np.abs(c))* np.exp( (1-int_const_1)/c ) * np.sqrt(2/np.pi) * (-int_const_2 + zeta) )), np.imag(( np.sqrt(np.abs(c))* np.exp( (1-int_const_1)/c ) * np.sqrt(2/np.pi) * (-int_const_2 + zeta) )) ))
        #print('input in with decimal erfinv: %.18f' %np.real( np.sqrt(np.abs(c))* np.exp( (1-int_const_1)/c ) * np.sqrt(2/np.pi) * (int_const_2 - zeta) ) )
        r = np.exp( (-1 + int_const_1 -c* scipy.special.erfinv(np.real( np.sqrt(np.abs(c))* np.exp( (1-int_const_1)/c ) * np.sqrt(2/np.pi) * (-int_const_2 + zeta) ))**2)/c )
    return r
    
def calc_r3_vp(zeta, c, int_const_1, int_const_2): #here I assumed v = - blabla because then you end up with this form
    
    if np.real(c)<0:
        # IERFI_APPROX = calc_ierfi( np.real( np.sqrt(np.abs(c))* np.exp( (1-int_const_1)/c ) * np.sqrt(2/np.pi) * (int_const_2 - zeta) ) )
        # #print('IERFI_APPROX:  %f ' %IERFI_APPROX)
        # r = np.exp( (-1 + int_const_1 +c*( IERFI_APPROX**2))/c )
        r = np.exp( (-1 + int_const_1 -c* ierf_approx(( np.sqrt(c)* np.exp( (1-int_const_1)/c ) * np.sqrt(2/np.pi) * (int_const_2 - zeta) ))**2)/c )
    else:
        #print('INPUT IN ERFINV : %f imag %f' %(np.real(( np.sqrt(np.abs(c))* np.exp( (1-int_const_1)/c ) * np.sqrt(2/np.pi) * (int_const_2 - zeta) )), np.imag(( np.sqrt(np.abs(c))* np.exp( (1-int_const_1)/c ) * np.sqrt(2/np.pi) * (int_const_2 - zeta) )) ))
        #print('input in with decimal erfinv: %f' %np.real( np.sqrt(np.abs(c))* np.exp( (1-int_const_1)/c ) * np.sqrt(2/np.pi) * (int_const_2 - zeta) ) )
        r = np.exp( (-1 + int_const_1 -c* scipy.special.erfinv(np.real( np.sqrt(np.abs(c))* np.exp( (1-int_const_1)/c ) * np.sqrt(2/np.pi) * (int_const_2 - zeta) ))**2)/c )
    return r
    #scipy.special.erfinv
    #ierf_approx


def calc_int_const_1(r, c, v ):
    int_const_1 = c*r**2/2 + 1/np.sqrt(1-v**2)

    return int_const_1

def calc_int_const_1_withlog(r, c, v ):
    int_const_1 = c*np.log(r) + 1/np.sqrt(1-v**2)

    return int_const_1

def calc_int_const_2(zeta, r, c, int_const_1):
    int_const_2 = zeta + 1j * np.log( - (1j * np.sqrt(3*c) * r + \
    np.sqrt( -6*int_const_1 -3 * c * r**2 +np.log(4) ) ) / \
    (12*int_const_1 -4* np.log(2) ) ) / np.sqrt(c)

    return int_const_2

def calc_int_const_2_vp(zeta, r, c, int_const_1):
    int_const_2 = zeta - 1j * np.log( -( -1j *np.sqrt(c)*r + np.sqrt(-2 +2*int_const_1 - c*r**2 ) )/( 4*(-1+int_const_1 ) ) ) \
    / np.sqrt(c)
    return int_const_2

def calc_int_const_2_vm(zeta, r, c, int_const_1):
    int_const_2 = zeta - 1j * np.log( ( 1j *np.sqrt(c)*r + np.sqrt(-2 +2*int_const_1 - c*r**2 ) )/( 4*(-1+int_const_1 ) ) ) \
    / np.sqrt(c)
    return int_const_2

def calc_int_const_2_withlog_m(zeta, r, c, int_const_1):
    
    if np.real(c)<0: 
        int_const_2 = zeta - ( np.exp( (-1+int_const_1)/c ) *np.sqrt(np.pi/2)* scipy.special.erfi( ( np.sqrt( (-1 + int_const_1 - c*np.log(r) )/np.abs(c) ) ) ) )/np.sqrt(np.abs(c))  
    else:
        int_const_2 = zeta - ( np.exp( (-1+int_const_1)/c ) *np.sqrt(np.pi/2)* scipy.special.erf( ( np.sqrt( (-1 + int_const_1 - c*np.log(r) )/np.abs(c) ) ) ) )/np.sqrt(np.abs(c))  
    return int_const_2

def calc_int_const_2_withlog_p(zeta, r, c, int_const_1):
    
    if np.real(c)<0:
        int_const_2 = zeta + ( np.exp( (-1+int_const_1)/c ) *np.sqrt(np.pi/2)* scipy.special.erfi( ( np.sqrt( (-1 + int_const_1 - c*np.log(r) )/abs(c) ) ) ) )/np.sqrt(abs(c)) 
    else:
        int_const_2 = zeta + ( np.exp( (-1+int_const_1)/c ) *np.sqrt(np.pi/2)* scipy.special.erf( ( np.sqrt( (-1 + int_const_1 - c*np.log(r) )/abs(c) ) ) ) )/np.sqrt(abs(c)) 
    
    return int_const_2


def calc_v2_p(r, c, int_const_1):
    v = np.sqrt( -4 + (-2*int_const_1 + c*r**2)**2 )/( 2*int_const_1 - c*r**2 )  #plusminus
    
    return v


def calc_v3_p(r, c, int_const_1):
    v = np.sqrt( ( -1 +int_const_1 -c*np.log(r) )*( 1 +int_const_1 -c*np.log(r) ) )/(int_const_1 - c * np.log(r)) #plusminus
    
    return v

def calc_v3_m(r, c, int_const_1):
    v = -np.sqrt( ( -1 +int_const_1 -c*np.log(r) )*( 1 +int_const_1 -c*np.log(r) ) )/(int_const_1 - c * np.log(r)) 
    
    return v

def plot_hipace_Ex(zeta_pos):
    ExmBy_path = './DATA/field_ExmBy_000000.h5'
    By_path = './DATA/field_By_000000.h5'
    
    ExmBy_g3d2 = Grid3d(ExmBy_path)
    By_g3d2 = Grid3d(By_path)
    
    ExmBy = np.transpose(ExmBy_g3d2.read(x2=0.0))
    By = np.transpose(By_g3d2.read(x2=0.0))  

    Ex = ExmBy + By
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    idx =np.abs(By_g3d2.get_zeta_arr() - zeta_pos).argmin()
    
    return By_g3d2.get_x_arr(2), Ex[:,idx]

def calc_transversal_probability_density(r_max, nx, nb, ni, rbunch, lbunch):
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
        if(r_array[i] <= rbunch):
            E[i] = np.real(c1prime) * r_array[i] 
        else:
            E[i] = np.real(c3prime) / r_array[i] 
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.plot(np.append(-r_array[::-1], r_array), np.append(-(E[::-1]) ,E[:]))
    x, Ex =plot_hipace_Ex(-0.2) # input = zeta pos 
    ax.plot(x, Ex)
    
    h5lp = H5Plot()
    h5lp.inherit_matplotlib_line_plots(ax)
    h5lp.write(savepath + '/' + savename + '.h5')
    plt.show()
    #use zeta array for shaping the elctric field to get the same spacing etc.
    #print(Earray)
    
    ionization_rate =np.zeros(shape=np.shape(r_array))
    ionization_rate[:] = vcalc_ion_rate(E*E_0, 1, 13.659843449,13.659843449, 0,0 ) # complete formulae for hydrogen
    #print(ionization_rate[1])
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #print('Deltat is %0.3e' %deltat)
    #computing the ionization probability
    

    savename = 'ion_probability_test'
    ion_probability = np.zeros(shape=(len(bunch_array) -1, len(ionization_rate)))
    for i in range(1,2): #len(bunch_array)):
        ion_probability[i-1,:] = 1.0 - np.exp(-ionization_rate[:] * deltat[i] )
        ax.plot(np.append(-r_array[::-1], r_array), np.append((ion_probability[i-1,::-1]) ,ion_probability[i-1,:]), label=('prob at lb: ' + str(np.around(bunch_array[i], decimals = 1) ) ) )
        #ax.legend()
        
        
    h5lp = H5Plot()
    h5lp.inherit_matplotlib_line_plots(ax)
    h5lp.write(savepath + '/' + savename + '.h5')
    plt.show()



def calc_ionization_length(start_r_array, c, rbunch, zeta_array):
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
    deltat = np.abs(zeta_array[1] -zeta_array[0]) / omega_p
    
    vcalc_ion_rate = np.vectorize(calc_ion_rate)
    
    if(start_r_array <= rbunch):
        E = np.real(c) * start_r_array * E_0
    else:
        E = np.real(c) / start_r_array * E_0
    
    #use zeta array for shaping the elctric field to get the same spacing etc.
    Earray = np.zeros(shape=np.shape(zeta_array))
    Earray[:] = E
    #print(Earray)
    
    ionization_rate =np.zeros(shape=np.shape(zeta_array))
    ionization_rate[:] = calc_ion_rate(Earray[1], 1, 13.659843449,13.659843449, 0,0 ) # complete formulae for hydrogen
    #print(ionization_rate[1])
    
    #print('Deltat is %0.3e' %deltat)
    #computing the ionization probability
    ion_probability = 1.0 - np.exp(-ionization_rate * deltat )
    
    #print(ion_probability[1])
    
    cum_prob = (1 - np.cumprod(1-ion_probability[::-1]))

    # print('cum prob')
    # print(cum_prob)

    k = -1
    for i in range(len(cum_prob)):
        if cum_prob[i]>0:
            k = i
            break
    
    return k
    
def calc_analytical_solution(start_r_array, nb, ni, lbunch, zeta_end, rbunch):
        shielding_factor = (1-start_r_array**2/rbunch**2)
        #forcing numpy to use complex squareroot
        c1 = 1/2*(shielding_factor*ni - nb) + 0j # 1/2*((1-start_r_array**2/rbunch**2) * ni - nb) + 0j ##
        c2 = 1/2*ni*shielding_factor + 0j # 1/2*((1-start_r_array**2/rbunch**2) *ni) + 0j # 
        c3 = 1/2*rbunch**2*( shielding_factor*ni - nb) + 0j# 1/2*rbunch**2*((1-start_r_array**2/rbunch**2) * ni - nb) + 0j
        c4 = 1/2*rbunch**2*( shielding_factor*ni) + 0j# 1/2*rbunch**2*((1-start_r_array**2/rbunch**2) * ni) + 0j
        #vectorizing the previous functions
        vcalc_r = np.vectorize(calc_r2_vm)
        vcalc_int_const_1 = np.vectorize(calc_int_const_1)
        vcalc_int_const_2 = np.vectorize(calc_int_const_2_vm)
        vcalc_v = np.vectorize(calc_v2_p) #used vm here
        vcalc_ion_length = np.vectorize(calc_ionization_length)
        
        #calculate integration constant 1 within the bunch at zeta = 0, v = 0:
        int_const_1 = vcalc_int_const_1(start_r_array, c1, 0)
        int_const_2 = vcalc_int_const_2(0, start_r_array, c1, int_const_1)
        #create zeta array:
        zeta_array = -np.arange(0, lbunch+0.01, 0.01)
        zeta_array_2 = -np.arange(lbunch, zeta_end, 0.01)
        
        
        
        
        if start_r_array<=rbunch:
            ion_length = calc_ionization_length(start_r_array, c1, rbunch, zeta_array)
        else:
            ion_length = calc_ionization_length(start_r_array, c3, rbunch, zeta_array)
        print('ion length:')
        print(ion_length)
        final_r_array = np.zeros(len(zeta_array))
    
        #start calculating the radius with along the zeta array
        if ion_length < 0:
            zeta_array = np.append(zeta_array, zeta_array_2)
            r_array2 = np.zeros(len(zeta_array_2))
            final_r_array = np.append(final_r_array,r_array2 )
            return zeta_array, np.real(final_r_array)
            
        for i in range(len(zeta_array)):
            if i<= ion_length:
                final_r_array[i] = 0 #start_r_array
            else:
            #Check, if we are still in the linear regime or if we have to switch to the 1/r potential
                if np.real(calc_r2_vm(zeta_array[i-ion_length], c1, int_const_1, (int_const_2))) < rbunch :
                    #print('linear field case')
                    final_r_array[i] = np.real(calc_r2_vm(zeta_array[i-ion_length], c1, int_const_1, (int_const_2)))
                else:
                    #print("1/r field case")
                    final_r_array[i] = np.real(calc_r2_vm(zeta_array[i-ion_length], c1, int_const_1, (int_const_2)))
                    if i == 0:
                        print('error particle starts outside rbunch!')
                        break
                    w = i
                    #calc v and the integration constants
                    v = calc_v2_p(final_r_array[w-1], c1, int_const_1)
                    int_const_1 = calc_int_const_1_withlog(final_r_array[w-1], c3, v )
                    int_const_2 = calc_int_const_2_withlog_p(zeta_array[w-1-ion_length], final_r_array[w-1], (c3), int_const_1 )
                    #Finish the rest of the bunch in the 1/r field regime
                    for j in range(w+1, len(zeta_array)):
                        final_r_array[j] = np.real(calc_r3_vm(zeta_array[j-ion_length], (c3), (int_const_1), (int_const_2) ) )
                    break
        
        # extend the r array for behind the beam
        r_array2 = np.zeros(len(zeta_array_2))
        final_r_array = np.append(final_r_array,r_array2 )
        
        #pick the right regime, depending on the starting position
        
        #case we are still in the linear regime
        if final_r_array[len(zeta_array)-1 ] < rbunch and final_r_array[len(zeta_array)-1 ]>0:
            #calc v
            v = calc_v2_p(final_r_array[len(zeta_array)-1 ], c1, int_const_1)
            # calc new int_const_1 with v at the end of the bunch:
            int_const_1 = calc_int_const_1(final_r_array[len(zeta_array)-1 ], c2, np.abs(v) )
            int_const_2 = calc_int_const_2_vm(-lbunch, final_r_array[len(zeta_array)-1 ], c2, int_const_1)
        
            for i in range(len(zeta_array_2)):
                
                # Again, check in which regime the particle is and adapt accordingly
                if np.real(calc_r2_vm(zeta_array_2[i], c2, int_const_1, (int_const_2))) < rbunch :
                    #print('linear field regime')
                    final_r_array[len(zeta_array)+ i] = np.real(calc_r2_vm(zeta_array_2[i], c2, int_const_1, (int_const_2)))
                    v = calc_v2_p(final_r_array[len(zeta_array)+ i], c2, int_const_1)
                
                else:
                    #print("1/r field regime")
                    w = i
                
                    #calc v and the integration constants
                    v = calc_v2_p(final_r_array[len(zeta_array)+ w-1], c2, int_const_1)
                    int_const_1 = calc_int_const_1_withlog(final_r_array[len(zeta_array)+ w-1], c4, v )
                    int_const_2 = calc_int_const_2_withlog_m(zeta_array_2[w-1], final_r_array[len(zeta_array)+ w-1], (c4), int_const_1 )
                
                    for j in range(w-1, len(zeta_array_2)):
                        if np.real(calc_r3_vm(zeta_array_2[j], (c4), (int_const_1), (int_const_2) ) ) > rbunch or j==0:
                            final_r_array[len(zeta_array)+ j] = np.real(calc_r3_vm(zeta_array_2[j], (c4), (int_const_1), (int_const_2) ) )
                        elif ((np.real(calc_r3_vm(zeta_array_2[j], (c4), (int_const_1), (int_const_2) ) ) < rbunch) and (np.real(calc_r3_vm(zeta_array_2[j-1], (c4), (int_const_1), (int_const_2) ) ) > rbunch and (j>=1 ))):
                            k = j
                            v = calc_v3_m(final_r_array[len(zeta_array) + k-1 ], c4, int_const_1)
                            int_const_1 = calc_int_const_1(final_r_array[len(zeta_array)+ k-1], c2, v )
                            int_const_2 = calc_int_const_2_vp(zeta_array_2[k-1], final_r_array[len(zeta_array)+ k-1], (c2), int_const_1 )
                        
                            for l in range(k-1, len(zeta_array_2)):
                                final_r_array[len(zeta_array) + l] = np.real(calc_r2_vm(zeta_array_2[l], (c2), (int_const_1), (int_const_2) ) )
                            break
                    break

        else: # <- case were we transition in the 1/r regime behind the bunch
            v = calc_v3_p(final_r_array[len(zeta_array)-1 ], c3, int_const_1)
            #print("v wert bei transition: %f " %v)
            int_const_1 = calc_int_const_1_withlog(final_r_array[len(zeta_array)-1 ], c4, np.abs(v) )
            int_const_2 = calc_int_const_2_withlog_m(-lbunch, final_r_array[len(zeta_array)-1 ], c4, int_const_1)
            #print('r wert bei transition:  %f' %(np.real(calc_r3_vm(zeta_array_2[0], (c4), (int_const_1), (int_const_2) ) )) )
            for i in range(len(zeta_array_2)):
                # check if the particle is coming back into the linear field regime 
                if np.real(calc_r3_vm(zeta_array_2[i], (c4), (int_const_1), (int_const_2) ) ) > rbunch: # or j==0:
                    final_r_array[len(zeta_array)+ i] = np.real(calc_r3_vm(zeta_array_2[i], (c4), (int_const_1), (int_const_2) ) )
                elif ((np.real(calc_r3_vm(zeta_array_2[i], (c4), (int_const_1), (int_const_2) ) ) < rbunch) and (np.real(calc_r3_vm(zeta_array_2[i-1], (c4), (int_const_1), (int_const_2) ) ) > rbunch and (i>=1 ))):
            
                    k = i
                    v = calc_v3_m(final_r_array[len(zeta_array) + k-1 ], c4, int_const_1)
                    int_const_1 = calc_int_const_1(final_r_array[len(zeta_array)+ k-1], c2, v )
                    int_const_2 = calc_int_const_2_vp(zeta_array_2[k-1], final_r_array[len(zeta_array)+ k-1], (c2), int_const_1 )
            
                    for l in range(k-1, len(zeta_array_2)):
                        final_r_array[len(zeta_array) + l] = np.real(calc_r2_vm(zeta_array_2[l], (c2), (int_const_1), (int_const_2) ) )
                    break
                else:
                    print('something went wrong in the case, where the particle was in the 1/r regime in and behind the beam and tried to come back to the linear regime')
            # 
        zeta_array = np.append(zeta_array, zeta_array_2)
        
        return zeta_array, np.real(final_r_array)



def main():
    
    modnum = 3
    getcontext().prec = 36
    basic_cols=['#75b765', '#808080', '#ffd700']
    basic_cols=['#2020ff' ,'#808080', '#ff2020']
    basic_cols=['#0000ff' , '#00ffff','#808080', '#ffff00', '#ff0000']
    my_cmap=LinearSegmentedColormap.from_list('mycmap', basic_cols)
    # display_cmap(my_cmap)
    
    parser = binSlab_parser() 
    args = parser.parse_args()

    # NHC=2
    proc_suffix_str = '_proc_'
    ppart_track_str = 'ppart_track'
    bin_fending = '.bin'
    
    
    zeta_min = -12
    #zeta_min = -12
    zeta_max = 4
    zeta_gridpoints = 1200
    # zeta_min = -8
    # zeta_max = 4
    # zeta_gridpoints = 300
    zeta_array = np.arange(zeta_min, zeta_max, abs(zeta_max - zeta_min)/zeta_gridpoints)
    numerical_data = np.load("/Users/diederse/desy/PIC-sim/HiPACE/tests/return_of_e2/rp_1_rb_1_2/data_numerical.npz")
    print(numerical_data.files)
    datan_zeta = numerical_data['arr_0']
    datan = numerical_data['arr_1']
    print(np.shape(datan_zeta))
    print(np.shape(datan))
    print(len(datan_zeta))
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
        
        # for k in range(len(w)):
        # 
        # 
        #     fig = plt.figure()
        #     if args.twodproj:
        #         ax = fig.add_subplot(111)
        #     else:
        #         ax = fig.add_subplot(111, projection='3d')
        #     ''' get min and max value for universal colorbar later '''
        #     if args.track_color == "u_tot":
        #         cmin = min(w[k][:,8])
        #         cmax = max(w[k][:,8])
        #         chosencmap = cm.jet
        #     elif args.track_color == "beta_z":
        #         cmin = -1
        #         cmax = 1
        #         chosencmap = my_cmap
        #     elif args.track_color == "beta_y":
        #         cmin = -1
        #         cmax = 1
        #         chosencmap = my_cmap
        #     else:
        #         print("This attribute doesn't exist or is not implemented yet")
        #         break
        #     # ##################### EXCLUDE HIPACE TRACKS FOR THE MOMENT
        #     ''' Splitting the array into each particle trajectory by proc tag'''
        #     d= np.split(w[k], np.where(np.diff(w[k][:, 7]) != 0)[0]+1) 
        # 
        #     for i in range(len(d)):
        # 
        #         ''' Splitting the array into each particle trajectory by part tag'''
        #         e = np.split(d[i], np.where(np.diff(d[i][:, 6]) != 0)[0]+1) 
        #         #print(np.shape(e))
        #         #print('particle tags %i (always look at 1 and ignore 0) %i'% (i, e[16][1, 6]))
        #         for j in range(int(np.floor(len(e)/modnum))):
        #             x=e[modnum*j][:,0]
        #             y=e[modnum*j][:,1]
        #             z=e[modnum*j][:,5] # IN CASE OF UNIONIZED PLASMA USE THIS TERM
        #             #starting_positions.append(x[ zeta_gridpoints -1]) #799])#
        #             z=zeta_array[zeta_gridpoints-len(y):] # IN CASE OF preionized PLASMA USE THIS TERM
        #             z=zeta_array[:len(y)]
        #             if args.track_color == "u_tot":
        #                 c=e[modnum*j][:,8]
        #             elif args.track_color == "beta_z":
        #                 c=e[modnum*j][:,4]/np.sqrt(1+e[modnum*j][:,8]**2)
        #             elif args.track_color == "beta_y":
        #                 c=e[modnum*j][:,2]/np.sqrt(1+e[modnum*j][:,8]**2)
        #             # print(len(x))
        #             # print(len(y))
        #             print("laenge track proc tag %i part tag %i ist %i" %(i, modnum*j, len(z)))
        # 
        #             # if args.twodproj:
        #             #     plot_2D_colourline(z,x,c, cmin, cmax)
        #             # else:
        #             #     plot_3D_colourline(z,y,x,c, cmin, cmax)
        # 
        #             # ################### TAKEN OUT TO FASTEN EVERYTHING FOR THE 1/R analysis!
        #             # if args.twodproj:
        #             #     if args.track_color == "u_tot":
        #             #         plot_2D_colourline(z,x,c, cmin, cmax)
        #             #     elif args.track_color == "beta_z":
        #             #         plot_2D_colourline_beta(z,x,c)
        #             #     elif args.track_color == "beta_y":
        #             #         plot_2D_colourline_beta(z,x,c)
        #             # else:
        #             #     if args.track_color == "u_tot":
        #             #         plot_3D_colourline(z,y,x,c, cmin, cmax)
        #             #     elif args.track_color == "beta_z":
        #             #         plot_3D_colourline_beta(z,y,x,c)
        #             #     elif args.track_color == "beta_y":
        #             #         plot_3D_colourline_beta(z,y,x,c)
        #             ###### ## ##ax.plot(z, y, x, label='particle track')
        # 
        #     if args.twodproj:
        # 
        # 
        # 
        #         input_r_array = np.array([ 0, 0.0941176488995552, 0.1882352977991104, 0.2823529541492462, 0.3764705955982208, 0.47058823704719543, 0.5647059082984924, 0.658823549747467, 0.7529411911964417, 0.8470588326454163, 0.9411764740943909, 1.0352941751480103, 1.1294118165969849, 1.2235294580459595, 1.317647099494934, 1.4117647409439087, 1.5058823823928833, 1.600000023841858, 1.6941176652908325, 1.7882353067398071, 1.8823529481887817, 1.9764705896377563])
        #         input_r_array = np.array([0.02346041053533554, 0.04692082107067108, 0.07038123160600662, 0.09384164214134216, 0.1173020526766777, 0.14076246321201324, 0.16422288119792938, 0.18768328428268433, 0.21114370226860046, 0.2346041053533554, 0.25806450843811035, 0.2815249264240265])
        #         #input_r_array = np.array([0.3128054738044739, 0.703812301158905, 1.094819188117981, 1.485826015472412, 1.8768328428268433, 2.2678396701812744, 2.658846616744995, 3.0498533248901367, 3.4408602714538574, 3.831866979598999, 4.222873687744141, 4.613880634307861]) #rp 5 mod 1
        #         input_r_array = np.array([ 0.03519061580300331, 0.07038123160600662, 0.10557185113430023, 0.14076246321201324, 0.17595307528972626, 0.21114370226860046, 0.24633431434631348]) #rp 0.3 mod 3
        #         #input_r_array = np.array([0.033,0.034,0.03516, 0.03517,0.03519061580300331]) #,0.0352])
        #         # fig = plt.figure()
        #         # if args.twodproj:
        #         #     ax = fig.add_subplot(111)
        #         # else:
        #         #     ax = fig.add_subplot(111, projection='3d')
        #         # 
        #         # for i in range(int(np.floor(len(datan_zeta)/modnum))):
        #         #     ax.plot(datan_zeta[modnum*i,:], datan[modnum*i,:], color = '#00cc00', linestyle = '--') #'#551a8b'
        #         # # 
        #         # # 
        #         modnum = 1
        #         #for i in range(int(np.floor(len(input_r_array)/modnum))):
        #         for i in range(len(input_r_array)):
        #             zeta_array3, r_matrix2 = calc_analytical_solution(input_r_array[i], 5, 1, 2, 12, 0.3) #36, 0.3)#8,2)  ###(start_r_array, nb, ni, lbunch, zeta_end, rbunch):
        #             ax.plot(zeta_array3[np.where(r_matrix2 >0)], r_matrix2[np.where(r_matrix2 >0)], color = 'black',linestyle = '-.' ) ##00cc00
        #         # zeta_array2, r_matrix = calc_analytical_solution( 1.999999999, 1.2, 1, 2, 8,2) #1.8823529481887817 1.9764705896377563
        #         # ax.plot(zeta_array2, r_matrix, 'r' )  
        #     else:
        #         ax = plt.gca()
        #         # ax.plot(data_zeta, np.zeros(len(data_zeta)), data2)
        # 
        #     ''' Set colorbar ''' 
        #     norm = matplotlib.colors.Normalize(
        #     vmin=np.min(cmin),
        #     vmax=np.max(cmax))
        # 
        #     # choose a colormap
        #     c_m = chosencmap
        # 
        #     # create a ScalarMappable and initialize a data structure
        #     s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
        #     s_m.set_array([])
        # 
        #     cbar = plt.colorbar(s_m)
        #     if args.track_color == "u_tot":
        #         cbar.ax.set_ylabel(r'$|u|$')
        #     elif args.track_color == "beta_z":
        #         cbar.ax.set_ylabel(r'$\beta_z$')
        #     elif args.track_color == "beta_y":
        #         cbar.ax.set_ylabel(r'$\beta_y$')
        # 
        # 
        # 
        #     #ax.set_xlim(-8, 0)
        #     # # ax.set_xlim(0, 300)
        #     #ax.set_ylim(-1/2, 6)
        #     ax.grid()
        #     if not args.twodproj:
        #         ax.set_zlabel(' x ')
        #         #ax.set_zlim(-6, 6)
        #     ax.set_xlabel(r'$\zeta$')
        #     ax.set_ylabel(' y ')
        # 
        #     numerical_solutions = mlines.Line2D([], [], color='#551a8b', label='numerical solutions')
        #     analytical_solutions = mlines.Line2D([], [], color='#00cc00', label='analytical solutions')
        #     #ax.legend(handles=[numerical_solutions,analytical_solutions ])
        # 
        #     print(starting_positions)
        #     plt.show()

    calc_transversal_probability_density(12, 256, 5, 1, 0.3, 2) #parameters from hipace.c ##(r_max, nx, nb, ni, rbunch, lbunch)


        #     # cax = plt.plot( array)
        #     # 
        #     # ax = plt.gca()
        #     # ax.set_ylabel(r'$\eta (\zeta)$', fontsize=14)
        #     # ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
        #     # 
        #     # if args.ifshow:
        #     #   plt.show()
        #     # 
        #     # head, savename = os.path.split(filepath)
        #     # 
        #     # mkdirs_if_nexist(args.savepath)
        #     # 
        #     # print('Saving: ' + args.savepath + '/' + savename + '.eps')
        #     # 
        #     # fig.savefig( args.savepath + '/' + savename + '.eps',
        #     #         format='eps')    

if __name__ == "__main__":
    main()  