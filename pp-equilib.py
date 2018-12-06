#!/usr/bin/env python3
import os
import sys
import math
import argparse
import numpy as np
# import importlib.util
import scipy
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import matplotlib
import matplotlib.pyplot as plt
# import matplotlib.colors as colors
# from matplotlib.colors import LogNorm
# from matplotlib import cm

pp_path = os.environ['PP_PATH']
sys.path.append(pp_path)
from pp_h5dat import Grid3d
from pp_h5dat import mkdirs_if_nexist
from pp_plt_tools import saveas_eps_pdf
from pp_plt_tools import saveas_png

def parseargs():

    desc='Text.'

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(  'nb_path',
                          metavar='NB_PATH',                          
                          type=str,
                          help = 'Path to nb file.')
    parser.add_argument(  '--wx-path',
                          dest='wx_path',
                          metavar='WX_PATH',
                          type=str,
                          help = 'Path to Wr file.')
    parser.add_argument(  'zeta_pos',
                          metavar='ZETA_POS',                          
                          type=float,
                          help = 'Zeta position.')      
    parser.add_argument(  '--mpsi-path',
                          dest='mpsi_path',        
                          metavar='MPSI_PATH',
                          type=str,
                          help = 'Path to mpsi file.')                            
    parser.add_argument(  '-v', '--verbose',
                          action='store_true',
                          dest='verbose',
                          default=True,
                          help = 'Print info (Default).')
    parser.add_argument(  '-q', '--quiet',
                          action='store_false',
                          dest='verbose',
                          help = 'Don''t print info.')

    return parser


class Data_xy_slice:
    def __init__(self, h5file, zeta_pos):
        self.__h5file = h5file
        self.__zeta_pos = zeta_pos

    def read(self):
        g3d = Grid3d(self.__h5file)
        self.__fxy = g3d.read(x0=self.__zeta_pos)
        self.__x_array = g3d.get_x_arr(1)
        self.__y_array = g3d.get_x_arr(2)

        return self.__x_array, self.__y_array, self.__fxy

    def plot_xy_slice(self):
        fig = plt.figure()
        plt.pcolor(self.__x_array, self.__y_array, self.__fxy)
        ax = plt.gca()        
        ax.set_xlabel(r'$k_p x$', fontsize=14)
        ax.set_ylabel(r'$k_p y$', fontsize=14)   

        plt.show()

    def cart_to_r_transform(self, x0 = 0.0, y0 = 0.0):
        if not (any(self.__x_array) or any(self.__y_array) or any(self.__fxy)):
            x, y, fxy = self.read()
        else:
            x, y, fxy = self.__x_array, self.__y_array, self.__fxy

        ix0 = np.argmin(np.abs(x-x0))
        if x[ix0] == x0:
            r = x[ix0:] - x0
        else:
            r = np.insert( x[(ix0+1):] - x0, 0, 0.0)

        f_integral = np.zeros(np.size(r),dtype=np.float32)

        f_interp = interpolate.interp2d(x,y,fxy)

        Nphi = np.size(x)
        phi_array = np.linspace(0,2*math.pi,num=Nphi)

        # Nphi = 1
        # phi_array = [0*math.pi/2 + 0.2]

        for phi in phi_array:
            x_interp = x0 + math.cos(phi) * r
            y_interp = y0 + math.sin(phi) * r

            if 0.0 <= phi < math.pi/2:
                f_integral += np.diagonal(f_interp(x_interp,y_interp))
            elif math.pi/2 <= phi < math.pi:
                f_integral += np.flip(np.diagonal(np.flipud(f_interp(x_interp,y_interp))))
            elif math.pi <= phi < 3*math.pi/2:
                f_integral += np.diagonal(np.flipud(np.fliplr(f_interp(x_interp,y_interp))))
            elif 3*math.pi/2 <= phi <= 2*math.pi:
                f_integral += np.flip(np.diagonal(np.fliplr(f_interp(x_interp,y_interp)))) 

        fr = f_integral/Nphi

        self.__r = r
        self.__fr = fr

        return r, fr

    def get_r_fr(self):
        return self.__r, self.__fr


# class n_xy_slice(Data_xy_slice):
    
# class Psi_xy_slice(Data_xy_slice):



class Equilib:
    def __init__(self, nr, Wr):
        pass


def main():
    parser = parseargs()
    args = parser.parse_args()

    nbdat = Data_xy_slice(args.nb_path, args.zeta_pos)
    nbdat.read()
    r_n, nbr = nbdat.cart_to_r_transform()

    psidat = Data_xy_slice(args.mpsi_path, args.zeta_pos)
    psidat.read()
    r_psi, mpsi = psidat.cart_to_r_transform()

    fig = plt.figure()
    plt.plot(r_n, nbr)
    ax = plt.gca()        
    ax.set_xlabel(r'$k_p r$', fontsize=14)
    ax.set_ylabel(r'$n_b$', fontsize=14)   

    fig = plt.figure()
    plt.plot(r_psi, mpsi)
    ax = plt.gca()        
    ax.set_xlabel(r'$k_p r$', fontsize=14)
    ax.set_ylabel(r'$W_x/E_0$', fontsize=14) 

    plt.show()

if __name__ == "__main__":
    main()
