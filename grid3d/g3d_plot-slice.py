#!/usr/bin/env python3

import numpy as np
from optparse import OptionParser
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import picdefs
from grid3d import Grid3d

def main():

  NUM_ARGS = 1

  usage = "usage: %prog [options] <file>"
  parser = OptionParser(usage=usage)
  (options, args) = parser.parse_args()

  if len(args) != NUM_ARGS:
    parser.error("This script requires an argument!")  

  # File 1
  g3d = Grid3d(picdefs.code.HIPACE)
  g3d.read(args[0])
  
  fig_zeta = plt.figure()
  plt.pcolormesh( g3d.get_zeta_arr(), 
                  g3d.get_x_arr(1), 
                  g3d.data[:,:,64].T, 
                  cmap='PuBu_r')
  ax = plt.gca()
  ax.set_ylabel(r'$k_p x$', fontsize=14)
  ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
  fig_zeta.savefig('plot_x_zeta.png', format='png')
  
  fig_xi = plt.figure()
  plt.pcolormesh( g3d.get_xi_arr(), 
                  g3d.get_x_arr(1), 
                  g3d.data[:,:,64].T, 
                  cmap='PuBu_r')
  ax = plt.gca()
  ax.set_ylabel(r'$k_p x$', fontsize=14)
  ax.set_xlabel(r'$k_p \xi$', fontsize=14)
  fig_xi.savefig('plot_x_xi.png', format='png')
  
if __name__ == "__main__":
    main()