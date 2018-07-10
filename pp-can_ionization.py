#!/usr/bin/env python3

import sys
import os
import math
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
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
  print('kp = ' + str(kp))



  #home path Path definitions
  Ez_path2 = '/Users/diederse/desy/PIC-sim/HiPACE/tests/the_return_of_the_local_electron/theo_test/DATA/field_Ez_000000.h5'
  ExmBy_path2 = '/Users/diederse/desy/PIC-sim/HiPACE/tests/the_return_of_the_local_electron/theo_test/DATA/field_ExmBy_000000.h5'
  EypBx_path2 = '/Users/diederse/desy/PIC-sim/HiPACE/tests/the_return_of_the_local_electron/theo_test/DATA/field_EypBx_000000.h5'
  Bx_path2 = '/Users/diederse/desy/PIC-sim/HiPACE/tests/the_return_of_the_local_electron/theo_test/DATA/field_Bx_000000.h5'
  By_path2 = '/Users/diederse/desy/PIC-sim/HiPACE/tests/the_return_of_the_local_electron/theo_test/DATA/field_By_000000.h5'
  

 # neutral_density_path2 = '/Users/diederse/desy/PIC-sim/HiPACE/tests/the_return_of_the_local_electron/theo_test/DATA/density_ionized_electrons_plasma_H_000000.h5'

  
  #Ez_path = '/Users/diederse/mountpoints/ed_scratch/simulations/hipace/tests/neutral1/DATA/field_Ez_000000.h5'
  #ExmBy_path = '/Users/diederse/mountpoints/ed_scratch/simulations/hipace/tests/neutral1/DATA/field_ExmBy_000000.h5'
  #EypBx_path = '/Users/diederse/mountpoints/ed_scratch/simulations/hipace/tests/neutral1/DATA/field_EypBx_000000.h5'
  #Bx_path = '/Users/diederse/mountpoints/ed_scratch/simulations/hipace/tests/neutral1/DATA/field_Bx_000000.h5'
  #By_path = '/Users/diederse/mountpoints/ed_scratch/simulations/hipace/tests/neutral1/DATA/field_By_000000.h5'
  
  #neutral_density_path = '/Users/diederse/mountpoints/ed_scratch/simulations/hipace/tests/neutral1/DATA/density_plasma_neutral_000000.h5'
  
  
  #vektorized function for calc ion rate
  vectorfunc = np.vectorize(calc_ion_rate)
  
  Ez_g3d = Grid3d(Ez_path2)
  ExmBy_g3d = Grid3d(ExmBy_path2)
  EypBx_g3d = Grid3d(EypBx_path2)
  Bx_g3d = Grid3d(Bx_path2)
  By_g3d = Grid3d(By_path2)

 # neutral_density_g3d = Grid3d(neutral_density_path2)
  
  
  
  #Reading in the data
  Ez = np.transpose(Ez_g3d.read(x2=0.0))
  ExmBy = np.transpose(ExmBy_g3d.read(x2=0.0))
  EypBx = np.transpose(EypBx_g3d.read(x2=0.0))
  Bx = np.transpose(Bx_g3d.read(x2=0.0))
  By = np.transpose(By_g3d.read(x2=0.0))  

  
  # neutral_density = np.transpose(neutral_density_g3d.read(x2=0.0))
  # plt.pcolormesh(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), neutral_density, cmap=cm.Blues_r) #
  # cb = plt.colorbar() 
  # plt.clim(0,-3)
  # #plt.imshow(neutral_density)
  # plt.show()
  
  
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

  if not os.path.exists('plots'):
    os.makedirs('plots')
  if not os.path.exists('plots/pp-ionization/'):
    os.makedirs('plots/pp-ionization/')

  ## computing the Magnitude of E:

  E_magnitude = np.sqrt((ExmBy + By)**2 + (EypBx - Bx)**2 + Ez**2 )* E_0
  


  # vmax1 = np.max(ExmBy + By)
  # vmin1 = np.min(ExmBy + By)
  # midpoint1 = 1 - vmax1/(vmax1 + np.abs(vmin1))
  # shifted_cmap1 = shiftedColorMap(cm.seismic, midpoint=midpoint1, name='shifted') 
  # 
  # ### Plotting Ex from HiPACE 
  # plt.contourf(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), ExmBy + By, 200, cmap=shifted_cmap1) # here cmap reds, since it is not divering
  # plt.ylabel(r'$k_p x$', fontsize =14)
  # plt.xlabel(r'$k_p \zeta$', fontsize =14)
  # cb = plt.colorbar()
  # cb.set_label(label = r'$ Ex$ HiPACE', fontsize = 14)
  # plt.savefig('plots/pp-ionization/E_x.png')
  # plt.show()
  
  
  # Plotting Ey from HiPACE
  # vmax = np.max(EypBx - Bx)
  # vmin = np.min(EypBx - Bx)
  # midpoint = 1 - vmax/(vmax + np.abs(vmin))
  # shifted_cmap = shiftedColorMap(cm.seismic, midpoint=midpoint, name='shifted')
  # 
  # ### Plotting the magnitude 
  # print('length By_g3d.get_x_arr(2) = ' + str(len(By_g3d.get_x_arr(2))))
  # plt.contourf(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), EypBx - Bx, 200, cmap=shifted_cmap) # here cmap reds, since it is not divering
  # plt.ylabel(r'$k_p x$', fontsize =14)
  # plt.xlabel(r'$k_p \zeta$', fontsize =14)
  # cb = plt.colorbar()
  # cb.set_label(label = r'$ Ey$', fontsize = 14)
  # plt.savefig('plots/pp-ionization/E_y.png')
  # plt.show()
  
  nb = 1.2
  ni = 1
  can_field = np.zeros(shape=(256,300))
  can_field2 = np.zeros(shape=(256,300))
  
  length_x_step = 24/256
  
  for i in range(0, 128):
      can_field[127-i, 150:200] = (63/64)**2 * nb * 0.5 / (length_x_step * (i + 0.5))
      can_field[128+i, 150:200] = -(63/64)**2 * nb * 0.5 / (length_x_step * (i + 0.5)) 
      # Rp needs to be adjusted here as well:
      can_field2[127-i, 0:200] = -(63/64)**2 * ni * 0.5 / (length_x_step * (i + 0.5))
      can_field2[128+i, 0:200] = +(63/64)**2 * ni * 0.5 / (length_x_step * (i + 0.5))
      
      #can_field2[127-i, 0:200] = -(255/64)**2 * ni * 0.5 / (0.09375 * (i + 0.5))
      #can_field2[128+i, 0:200] = +(255/64)**2 * ni * 0.5 / (0.09375 * (i + 0.5))
      
      ### FOR Rp = 2
      can_field[127-i, 150:200] = (129/64)**2 * nb * 0.5 / (length_x_step * (i + 0.5))
      can_field[128+i, 150:200] = -(129/64)**2 * nb * 0.5 / (length_x_step * (i + 0.5)) 
      # Rp needs to be adjusted here as well:
      can_field2[127-i, 0:200] = -(129/64)**2 * ni * 0.5 / (length_x_step * (i + 0.5))
      can_field2[128+i, 0:200] = +(129/64)**2 * ni * 0.5 / (length_x_step * (i + 0.5))

      ### for rp=rb=2
  for i in range(0,22):
      can_field[127-i, 150:200] = (length_x_step * (i + 0.5)) * nb * 0.5 
      can_field[128+i, 150:200] = -length_x_step * (i + 0.5) * nb  * 0.5
      can_field2[127-i, 0:200] = -(length_x_step * (i + 0.5)) * ni * 0.5 
      can_field2[128+i, 0:200] = +length_x_step * (i + 0.5) * ni  * 0.5 
  
  
  # for i in range(0,11):
  #     can_field[127-i, 150:200] = (length_x_step * (i + 0.5)) * nb * 0.5 
  #     can_field[128+i, 150:200] = -length_x_step * (i + 0.5) * nb  * 0.5 
  
  #for Rp = 255/64
  # for i in range(0, 43):
  #     can_field2[127-i, 0:200] = -(0.09375 * (i + 0.5)) * ni * 0.5 
  #     can_field2[128+i, 0:200] = +0.09375 * (i + 0.5) * ni  * 0.5 

  #for Rp = 63/64
  # for i in range(0, 11):
  #     can_field2[127-i, 0:200] = -(length_x_step * (i + 0.5)) * ni * 0.5 
  #     can_field2[128+i, 0:200] = +length_x_step * (i + 0.5) * ni  * 0.5 

  
  #for Rp = 33/64
  # for i in range(0, 6):
  #    can_field2[127-i, 0:200] = -(0.09375 * (i + 0.5)) * ni * 0.5 
  #    can_field2[128+i, 0:200] = +0.09375 * (i + 0.5) * ni  * 0.5 

  can_field += can_field2
  
  plt.plot(np.arange(-12,12, length_x_step), can_field[:, 175])
  plt.show()
  
  vmax = np.max(can_field)
  vmin = np.min(can_field)
  midpoint = 1 - vmax/(vmax + np.abs(vmin))
  shifted_cmap = shiftedColorMap(cm.seismic, midpoint=midpoint, name='shifted')
  plt.contourf(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), can_field, 200, cmap=shifted_cmap) # here cmap reds, since it is not divering
  plt.ylabel(r'$k_p x$', fontsize =14)
  plt.xlabel(r'$k_p \zeta$', fontsize =14)
  cb = plt.colorbar()
  cb.set_label(label = r'$ E_r/E_0$ analytical', fontsize = 14)
  plt.axvline(x = -2, color='black')
  plt.axhline(y = 129/64, color='black')
  plt.text(-7.5, 11, 'IV', fontsize=16)
  plt.text(-1, 11, 'III', fontsize=16)
  plt.text(-1, 0.5, 'I', fontsize=16)
  plt.text(-7.5, 0.5, 'II', fontsize=16)
  plt.xlim(-8, 0)
  plt.ylim(0,12)
  plt.show()
  
  
  
  ## Computing the ionization rate:
  #ionization_rate = np.nan_to_num(omega_alpha / (2.0 * np.pi) *4.0 * E_alpha/ E_magnitude * np.exp (-2.0/3.0 * E_alpha/ E_magnitude ) )# this is just correct for hydrogen!
  ionization_rate = vectorfunc(abs(can_field)*E_0, 1, 13.659843449,13.659843449, 0,0 ) # complete formulae for hydrogen
  #ionization_rate = vectorfunc(E_magnitude, 2, 54.4177650, 24.58738880, 0,0 ) 
  
  potential = np.zeros(shape=(256,300))
  #potential[0:128, :] = -np.cumsum(can_field[0:128, :], axis=0)
  potential[128:, :] = -np.cumsum(can_field[128:, :]*0.09375, axis=0)
  # for i in range(128):
  #     potential[128+i] = potential[127-i,:]  
  
  
  startradius = 63/63
  endradius = 1.57
  def calc_pot_at_x_within_bunch(x):
      return potential[int(128+np.floor(x/0.09375)),190]

  def calc_pot_at_x_behind_bunch(x):
      return -potential[int(128+np.floor(x/0.09375)),10]

  
  print('potential at 1.567:')
  print(calc_pot_at_x_within_bunch(endradius))

  print('potential at 63/64:')
  print(calc_pot_at_x_within_bunch(startradius))
  
  print('therefore the total kinetic energy the electron gained is: ')
  print( calc_pot_at_x_within_bunch(endradius) - calc_pot_at_x_within_bunch(startradius) +  calc_pot_at_x_behind_bunch(endradius) )
  
  plt.plot(By_g3d.get_x_arr(2)[128:], -potential[128:,10], label=r'$E_{pot}$ behind beam analytical')
  plt.plot(By_g3d.get_x_arr(2)[128:], potential[128:,190], label=r'$-E_{pot}$ within beam analytical')
  #plt.plot(By_g3d.get_x_arr(2)[128:], can_field[128:,10], label=r'$E_r/E_0$ analytical')
  plt.legend()
  plt.grid()
  plt.show()
  
  vmax = np.max(potential)
  vmin = np.min(potential)
  midpoint = 1 - vmax/(vmax + np.abs(vmin))
  shifted_cmap = shiftedColorMap(cm.seismic, midpoint=midpoint, name='shifted')
  plt.contourf(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), potential, 200, cmap=shifted_cmap) # here cmap reds, since it is not divering
  plt.ylabel(r'$k_p x$', fontsize =14)
  plt.xlabel(r'$k_p \zeta$', fontsize =14)
  cb = plt.colorbar()
  cb.set_label(label = r'$ V$ analytical', fontsize = 14)
  
  plt.show()
  
  ekin_at_r1 = -(potential- potential[ 133, 190])
  vmax = np.max(ekin_at_r1)
  vmin = np.min(ekin_at_r1)
  midpoint = 1 - vmax/(vmax + np.abs(vmin))
  shifted_cmap = shiftedColorMap(cm.seismic, midpoint=midpoint, name='shifted')
  plt.contourf(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), ekin_at_r1, 200, cmap=shifted_cmap) # here cmap reds, since it is not divering
  plt.ylabel(r'$k_p x$', fontsize =14)
  plt.xlabel(r'$k_p \zeta$', fontsize =14)
  cb = plt.colorbar()
  cb.set_label(label = r'$ epot @ r = 1 $ analytical', fontsize = 14)
  
  plt.show()



  ### Plotting the magnitude 
  plt.contourf(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), E_magnitude, 200, cmap=cm.Reds) # here cmap reds, since it is not divering
  plt.ylabel(r'$k_p x$', fontsize =14)
  plt.xlabel(r'$k_p \zeta$', fontsize =14)
  cb = plt.colorbar()
  cb.set_label(label = r'$ |E|$ HiPACE', fontsize = 14)
  plt.savefig('plots/pp-ionization/E_magnitude.png')
  plt.show()
  
  #print('The maximum ionization rate is %0.3e' % np.max(ionization_rate))
  ### Plotting the ionization rate:
  plt.contourf(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), can_field, 200, cmap=cm.seismic) # here cmap reds, since it is not divering
  plt.ylabel(r'$k_p x$', fontsize =14)
  plt.xlabel(r'$k_p \zeta$', fontsize =14)
  cb = plt.colorbar()
  cb.set_label(label = r'$ |E| analytical$', fontsize = 14)
  plt.savefig('plots/pp-ionization/can_field.png')
  plt.show()
  
  plt.plot(By_g3d.get_x_arr(2), (ExmBy + By)[:,190], label=r'$Ex/E_0$ HiPACE')
  plt.plot(By_g3d.get_x_arr(2), can_field[:,190], label=r'$Ex/E_0$ analytical')
  plt.legend()
  plt.show()

  print('E at roughly 0.5 within beam: %f ' %can_field[122, 190])
  print('E at roughly 0.5 behind beam : %f ' %can_field[122, 1])
  print('E at roughly 4 within beam: %f ' %can_field[85, 190])
  print('E at roughly 4 behind beam : %f ' %can_field[85, 1])
  print('E at roughly 4 within beam: %f ' %can_field[85, 190])
  print('E at roughly 4 behind beam : %f ' %can_field[85, 1])
  #print('The maximum ionization rate is %0.3e' % np.max(ionization_rate))
  ### Plotting the ionization rate:
  plt.contourf(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), ionization_rate, 200, cmap=cm.Reds) # here cmap reds, since it is not divering
  plt.ylabel(r'$k_p x$', fontsize =14)
  plt.xlabel(r'$k_p \zeta$', fontsize =14)
  cb = plt.colorbar()
  cb.set_label(label = 'Ionization rate [1/s]', fontsize = 14)
  plt.savefig('plots/pp-ionization/ion_rate.png')
  plt.show()
  



  cbarvektor = np.linspace(0,1,11, endpoint = True)
  #print(By_g3d.get_zeta_arr())
  #print(np.abs(By_g3d.get_zeta_arr()[1] -By_g3d.get_zeta_arr()[0]))
  deltat = np.abs(By_g3d.get_zeta_arr()[1] -By_g3d.get_zeta_arr()[0]) / omega_p
  print('Deltat is %0.3e' %deltat)
  #computing the ionization probability
  ion_probability = 1.0 - np.exp(-ionization_rate * deltat )
  #print('The maximum ionization probability is %0.3e' % np.max(ion_probability))
  ## Plotting the ionization probability:
  plt.contourf(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), ion_probability, 200,cmap=cm.Reds) # here cmap reds, since it is not divering  vmin=0, vmax=1,
  plt.ylabel(r'$k_p x$', fontsize =14)
  plt.xlabel(r'$k_p \zeta$', fontsize =14)
  
  cb.set_label(label = 'Ionization probability', fontsize = 14)
  plt.show()

  # cbarvektor = np.linspace(0,1,11, endpoint = True)
  # #print(By_g3d.get_zeta_arr())
  # #print(np.abs(By_g3d.get_zeta_arr()[1] -By_g3d.get_zeta_arr()[0]))
  # deltat = np.abs(By_g3d.get_zeta_arr()[1] -By_g3d.get_zeta_arr()[0]) / omega_p
  # print('Deltat is %0.3e' %deltat)
  # #computing the ionization probability
  # ion_probability = 1.0 - np.exp(-ionization_rate * deltat )
  # #print('The maximum ionization probability is %0.3e' % np.max(ion_probability))
  # ### Plotting the ionization probability:
  # plt.contourf(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), ion_probability, 200,cmap=cm.Reds) # here cmap reds, since it is not divering  vmin=0, vmax=1,
  # plt.ylabel(r'$k_p x$', fontsize =14)
  # plt.xlabel(r'$k_p \zeta$', fontsize =14)
  # 
  # cb = plt.colorbar(ticks = cbarvektor)
  # #cb.set_clim(vmin=0,vmax=1)
  # #cb.set_ticks(cbarvektor)
  # #plt.clim(0,1)
  # cb.set_label(label = 'Ionization probability', fontsize = 14)
  # plt.show()

  '''
  #calculate the cumulative probability:
  cum_prob = (1 - np.cumprod(1-ion_probability[:,::-1], axis=1))[:,::-1] #np.flip(np.flip(ion_probability, 1).cumsum(), 1) #
  
  
  ### Plotting the ionization probability:
  plt.contourf(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), cum_prob, 200, cmap=cm.Reds) # here cmap reds, since it is not divering
  plt.ylabel(r'$k_p x$', fontsize =14)
  plt.xlabel(r'$k_p \zeta$', fontsize =14)
  cb = plt.colorbar(ticks = cbarvektor) 
  #plt.clim(0,1.01)
  cb.set_label(label = 'Ionization probability', fontsize = 14)
  plt.savefig('plots/pp-ionization/cum_ion_prob.png')
  plt.show()

  
  ### Plotting the ionization probability:
  plt.pcolormesh(By_g3d.get_zeta_arr(), By_g3d.get_x_arr(2), ion_probability, cmap=cm.Reds) #here cmap reds, since it is not divering
  plt.ylabel(r'$k_p x$', fontsize =14)
  plt.xlabel(r'$k_p \zeta$', fontsize =14)
  cb = plt.colorbar(ticks = cbarvektor) 
  plt.clim(0,1)
  cb.set_label(label = 'Ionization probability', fontsize = 14)
  plt.savefig('plots/pp-ionization/ion_prob.png')
  plt.show()

  
  
  

  
  #Plotting the differences between HiPACE and analytical model 
  differenz = ( cum_prob + 1 ) - neutral_density # cum prob + 1 because it is preionized to level 1
  
  
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
  plt.savefig('plots/pp-ionization/diff_ion_prob.png')
  plt.show()

  '''
  '''
  Idee, aus E_magnitude kann die Rate berechnet werden, aus der Rate die Wahrscheinlichkeit.
  Ueber eine kumultative Summe von Rechts nach links kann dann das Bild produziert werden, 
  dass mit dem Ionization modul entstehen sollte. 
  '''


if __name__ == "__main__":
    main() 

