#!/usr/bin/env python3

import sys
import os
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib.colors import LogNorm
import scipy
from pp_h5dat import Grid3d
from scipy import special


from mpl_toolkits.axes_grid1 import AxesGrid

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

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
    #print('n_eff = %f' % n_eff)
    #print('C_2 = %f' % C_2)
    #print('Quant_factor = %f' % quant_prefactor)
    value1 = 0;
    value2 = 0;
    value3 = 0;
    if (elec_field/E_alpha > 1e-3):
    
        value1 = omega_alpha * quant_prefactor * C_2 * ion_pot/ (2* HYDROGEN_IONIZATION_ENERGIE_EV) 
        value2 = (2*E_alpha/elec_field * (ion_pot/ HYDROGEN_IONIZATION_ENERGIE_EV)**1.5)**(2*n_eff -1) 
        value3 = np.exp(-2.0/3.0 * E_alpha/elec_field * (ion_pot/ HYDROGEN_IONIZATION_ENERGIE_EV)**1.5 )
    
    return value1*value2*value3 









def main():

  #### Defining constants as in HiPACE
  #define M_ELEC_PER_M_NUCL         5.4857990946e-4
  #define E_MUON_MASS_RATIO         4.83633170e-3 /* from NIST database */
  ELECTRON_CHARGE_IN_COUL   =    1.60217657e-19
  ELECTRON_MASS_IN_KG       =    9.10938291e-31
  SPEED_OF_LIGHT_IN_M_PER_S =    299792458
  VAC_PERMIT_FARAD_PER_M    =    8.854187817e-12
  HYDROGEN_IONIZATION_ENERGIE_EV = 13.659843449
  
  omega_alpha               =    4.13e16; # [s^-1]
  E_alpha                   =    5.1e11; # [V/m]
  elec_density              =    1e+23 # [1/m^3]
  omega_p = np.sqrt(4* np.pi * elec_density * (ELECTRON_CHARGE_IN_COUL**2)/ (VAC_PERMIT_FARAD_PER_M * ELECTRON_MASS_IN_KG)); # calculation of the plasma frequency
  #print('omega_p is %0.3e' %omega_p)
  E_0 = omega_p * ELECTRON_MASS_IN_KG * SPEED_OF_LIGHT_IN_M_PER_S / ELECTRON_CHARGE_IN_COUL;                        #calculation of the denormalization factor for the electric field
  #print('E_0 is %0.3e' %E_0)
  kp = SPEED_OF_LIGHT_IN_M_PER_S/omega_p
  #home path Path definitions
  Ez_path = './DATA/field_Ez_000000.0.h5'
  ExmBy_path = './DATA/field_ExmBy_000000.0.h5'
  EypBx_path = './DATA/field_EypBx_000000.0.h5'
  Bx_path = './DATA/field_Bx_000000.0.h5'
  By_path = './DATA/field_By_000000.0.h5'
  
  #neutral_density_path = '/Users/diederse/desy/PIC-sim/HiPACE/tests/interaction_tests2/Cs/DATA/density_plasma_neutral_000000.h5'
  neutral_density_path = './DATA/density_plasma_H_000000.0.h5'
  ### nersc path HAS TO BE MOUNTED
  
  #Ez_path = '/Users/diederse/mountpoints/ed_scratch/simulations/hipace/tests/neutral1/DATA/field_Ez_000000.h5'
  #ExmBy_path = '/Users/diederse/mountpoints/ed_scratch/simulations/hipace/tests/neutral1/DATA/field_ExmBy_000000.h5'
  #EypBx_path = '/Users/diederse/mountpoints/ed_scratch/simulations/hipace/tests/neutral1/DATA/field_EypBx_000000.h5'
  #Bx_path = '/Users/diederse/mountpoints/ed_scratch/simulations/hipace/tests/neutral1/DATA/field_Bx_000000.h5'
  #By_path = '/Users/diederse/mountpoints/ed_scratch/simulations/hipace/tests/neutral1/DATA/field_By_000000.h5'
  
  #neutral_density_path = '/Users/diederse/mountpoints/ed_scratch/simulations/hipace/tests/neutral1/DATA/density_plasma_neutral_000000.h5'
  
  
  #vektorized function for calc ion rate
  vectorfunc = np.vectorize(calc_ion_rate)
  
  Ez_g3d = Grid3d(Ez_path)
  ExmBy_g3d = Grid3d(ExmBy_path)
  EypBx_g3d = Grid3d(EypBx_path)
  Bx_g3d = Grid3d(Bx_path)
  By_g3d = Grid3d(By_path)

  neutral_density_g3d = Grid3d(neutral_density_path)
  
  
  
  #Reading in the data
  Ez = np.transpose(Ez_g3d.read(x2=0.0))
  ExmBy = np.transpose(ExmBy_g3d.read(x2=0.0))
  EypBx = np.transpose(EypBx_g3d.read(x2=0.0))
  Bx = np.transpose(Bx_g3d.read(x2=0.0))
  By = np.transpose(By_g3d.read(x2=0.0))  
  
  
  neutral_density = np.transpose(neutral_density_g3d.read(x2=0.0))
  plt.pcolormesh(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), neutral_density, cmap=cm.Blues) #
  cb = plt.colorbar() 
  #plt.clim(0,3)
  #plt.imshow(neutral_density)
  plt.show()
  plt.close(fig)
  
  '''
  #### Test to plot the By field 
  print(By_g3d.get_x_arr(0)) #'zeta array limits are %f' % 
  plt.contourf(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), By, 200, cmap=cm.seismic)
  plt.ylabel(r'$k_p x$', fontsize =14)
  #plt.xlim(By_g3d.get_zeta_arr)
  #plt.ylim(By_g3d.get_x_arr)
  plt.xlabel(r'$k_p \zeta$', fontsize =14)
  cb = plt.colorbar()
  cb.set_label(label = r'$B_y/E_0$', fontsize = 14)
  plt.show()
  '''

#check for directory 

  if not os.path.exists('./plots'):
    os.makedirs('./plots')
  if not os.path.exists('./plots/pp-ionization/'):
    os.makedirs('./plots/pp-ionization/')

  ## computing the Magnitude of E:
 
  E_magnitude = np.sqrt((ExmBy + By)**2 + (EypBx - Bx)**2 + Ez**2 )* E_0
  
  number_above_5_6_e10 =  np.where( E_magnitude > 5.7e10)
  number = E_magnitude[number_above_5_6_e10]
  print(number)
  print('the number of grid points with an electric field above 5.6e10 is %0.3e' % len(number))
  #print(number_above_5_6_e10)
  print('The maximum magnitude of the electric field is %0.3e' % np.max(E_magnitude))
  ### Plotting the magnitude 
  plt.contourf(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), E_magnitude, 200, cmap=cm.Reds) # here cmap reds, since it is not divering
  plt.ylabel(r'$k_p x$', fontsize =14)
  plt.xlabel(r'$k_p \zeta$', fontsize =14)
  cb = plt.colorbar()
  cb.set_label(label = r'$ |E|$', fontsize = 14)
  plt.savefig('./plots/pp-ionization/E_magnitude.png')
  plt.show()
  plt.close(fig)
  
  ## Computing the ionization rate:
  ionization_rate = np.nan_to_num(omega_alpha / (2.0 * np.pi) *4.0 * E_alpha/ E_magnitude * np.exp (-2.0/3.0 * E_alpha/ E_magnitude ) )# this is just correct for hydrogen!
  #ionization_rate = vectorfunc(E_magnitude, 1, 13.659843449,13.659843449, 0,0 ) # complete formulae for hydrogen
  #ionization_rate = vectorfunc(E_magnitude, 2, 54.4177650, 24.58738880, 0,0 ) 
  
  #print('The maximum ionization rate is %0.3e' % np.max(ionization_rate))
  ### Plotting the ionization rate:
  plt.contourf(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), ionization_rate, 200, cmap=cm.Reds) # here cmap reds, since it is not divering
  plt.ylabel(r'$k_p x$', fontsize =14)
  plt.xlabel(r'$k_p \zeta$', fontsize =14)
  cb = plt.colorbar()
  cb.set_label(label = 'Ionization rate [1/s]', fontsize = 14)
  plt.savefig('./plots/pp-ionization/ion_rate.png')
  plt.show()
  plt.close(fig)
  
  #plotting the rate vs electric field in general
  #E_magnitude2 = 10**np.array(list(range(10, 15)))
  #ionization_rate_plot = np.nan_to_num(omega_alpha / (2.0 * np.pi) *4.0 * E_alpha/ E_magnitude2 * np.exp (-2.0/3.0 * E_alpha/ E_magnitude2 ) )
  #plt.plot(E_magnitude2, ionization_rate_plot)
  #plt.show()

 
  cbarvektor = np.linspace(0,1,11, endpoint = True)
  #print(By_g3d.get_zeta_arr())
  #print(np.abs(By_g3d.get_zeta_arr()[1] -By_g3d.get_zeta_arr()[0]))
  deltat = np.abs(By_g3d.get_zeta_arr()[1] -By_g3d.get_zeta_arr()[0]) / omega_p
  print('Deltat is %0.3e' %deltat)
  #computing the ionization probability
  ion_probability = 1.0 - np.exp(-ionization_rate * deltat )
  #print('The maximum ionization probability is %0.3e' % np.max(ion_probability))
  ### Plotting the ionization probability:
  plt.contourf(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), ion_probability, 200,cmap=cm.Reds) # here cmap reds, since it is not divering  vmin=0, vmax=1,
  plt.ylabel(r'$k_p x$', fontsize =14)
  plt.xlabel(r'$k_p \zeta$', fontsize =14)

  cb = plt.colorbar(ticks = cbarvektor)
  #cb.set_clim(vmin=0,vmax=1)
  #cb.set_ticks(cbarvektor)
  #plt.clim(0,1)
  cb.set_label(label = 'Ionization probability', fontsize = 14)
  plt.show()
  plt.close(fig)

  #calculate the cumulative probability:
  cum_prob = (1 - np.cumprod(1-ion_probability[:,::-1], axis=1))[:,::-1] #np.flip(np.flip(ion_probability, 1).cumsum(), 1) #
  
  
  ### Plotting the ionization probability:
  plt.contourf(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), cum_prob, 200, cmap=cm.Reds) # here cmap reds, since it is not divering
  plt.ylabel(r'$k_p x$', fontsize =14)
  plt.xlabel(r'$k_p \zeta$', fontsize =14)
  cb = plt.colorbar(ticks = cbarvektor) 
  #plt.clim(0,1.01)
  cb.set_label(label = 'Ionization probability', fontsize = 14)
  plt.savefig('./plots/pp-ionization/cum_ion_prob.png')
  plt.show()
  plt.close(fig)
  
  ### Plotting the ionization probability:
  plt.pcolormesh(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), ion_probability, cmap=cm.Reds) #here cmap reds, since it is not divering
  plt.ylabel(r'$k_p x$', fontsize =14)
  plt.xlabel(r'$k_p \zeta$', fontsize =14)
  cb = plt.colorbar(ticks = cbarvektor) 
  plt.clim(0,1)
  cb.set_label(label = 'Ionization probability', fontsize = 14)
  plt.savefig('./plots/pp-ionization/ion_prob.png')
  plt.show()
  plt.close(fig)
  
  
  

  
  #Plotting the differences between HiPACE and analytical model 
  differenz = ( cum_prob  ) - neutral_density # cum prob + 1 because it is preionized to level 1
  
  
  #producing a asymmetric Colorbar
  vmax = np.max(differenz)
  vmin = np.min(differenz)
  midpoint = 1 - vmax/(vmax + np.abs(vmin))
  shifted_cmap = shiftedColorMap(cm.seismic, midpoint=midpoint, name='shifted')
  
  
  plt.contourf(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), differenz, 200, cmap=shifted_cmap ) # cm.seismic) # here cmap reds, since it is not diverging
  plt.ylabel(r'$k_p x$', fontsize =14)
  plt.xlabel(r'$k_p \zeta$', fontsize =14)
  cb = plt.colorbar() 
  plt.clim(vmin,vmax)
  cb.set_label(label = r'$N_{H^+ analyt} - N_{H^+ HiPace}$', fontsize = 14)
  plt.savefig('./plots/pp-ionization/diff_ion_prob.png')
  plt.show()
  plt.close(fig)
  print('summed absolut difference: %f' %(np.sum(abs(differenz)))) 
  '''
  Idee, aus E_magnitude kann die Rate berechnet werden, aus der Rate die Wahrscheinlichkeit.
  Ueber eine kumultative Summe von Rechts nach links kann dann das Bild produziert werden, 
  dass mit dem Ionization modul entstehen sollte. 
  '''


if __name__ == "__main__":
    main() 

