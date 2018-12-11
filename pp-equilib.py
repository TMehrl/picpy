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
from scipy.optimize import curve_fit
from scipy.ndimage import filters

import matplotlib
import matplotlib.pyplot as plt

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
    parser.add_argument(  'wx_path',
                          metavar='WX_PATH',
                          type=str,
                          default=None,                           
                          help = 'Path to Wr file.')
    parser.add_argument(  'zeta_pos', 
                          metavar='ZETA_POS',                          
                          type=float,
                          nargs='+',
                          default=None,                           
                          help = 'Zeta position.')      
    parser.add_argument(  '--mpsi-path',
                          dest='mpsi_path',        
                          metavar='MPSI_PATH',
                          type=str,
                          default=None,                          
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
    parser.add_argument(  '--rlim',
                          help='Upper limit of radius to be taken into account.',
                          action='store',
                          dest="rlim",
                          metavar='RLIM',
                          type=float,
                          default=None)
    parser.add_argument(  '--n-zeta-avg',
                          help='Numer of grid cells to longtidudinally average over.',
                          action='store',
                          dest="nzetaavg",
                          metavar='NAVG',
                          type=int,
                          default=None)    
    parser.add_argument(  '--ssmooth',
                          help='Number of gridpoints (sigma) to be Gaussian smoothed.',
                          action='store',
                          dest="smooth_sigma",
                          metavar='SIGMA',
                          type=int,
                          default=0)    
    parser.add_argument(  "-s", "--save-path",
                          action="store",
                          dest="savepath",
                          metavar="PATH",
                          default='./plots/equilib',
                          help = """Path to which generated files will be saved.
                              (default: %(default)s)""")
    return parser


class Data_xy_slice:
    def __init__(self, h5file, zeta_pos):
        self.__h5file = h5file
        self.__zeta_pos = zeta_pos

    def read(self, navg=None):
        g3d = Grid3d(self.__h5file)
        self.__fxy = g3d.read(x0=self.__zeta_pos, navg=navg)
        self.__x_array = g3d.get_x_arr(1)
        self.__y_array = g3d.get_x_arr(2)

        return self.__x_array, self.__y_array, self.__fxy


def exp_func(x, a, b):
    return a * np.exp(-b * x)

def cart_int(x,y,fxy,axis='x'):
    if axis == 'x':
        f_int = integrate.cumtrapz(fxy, x, axis=0, initial=0)
    elif axis == 'y':
        f_int = integrate.cumtrapz(fxy, y, axis=1, initial=0)
    else:
        print('ERROR: Only "x" or "y" allowed for axis!')
        sys.exit(1)
    return x, y, f_int

def cart_grad(x,y,fxy,axis='x'):
    if axis == 'x':
        f_grad = np.gradient(fxy, x[1]-x[0], axis=0)
    elif axis == 'y':
        f_grad = np.gradient(fxy, y[1]-y[0], axis=1)
    else:
        print('ERROR: Only "x" or "y" allowed for axis!')
        sys.exit(1)
    return x, y, f_grad

def plot_xy_slice(x,y,fxy):
    fig = plt.figure()
    plt.pcolor(x, y, fxy)
    ax = plt.gca()        
    ax.set_xlabel(r'$k_p x$', fontsize=14)
    ax.set_ylabel(r'$k_p y$', fontsize=14)   
    plt.colorbar()
    plt.show()

def cart_to_r_transform(x, y, fxy, x0 = 0.0, y0 = 0.0, rlim = None):

    ix0 = np.argmin(np.abs(x-x0))

    if rlim != None:
        ix1 = np.argmin(np.abs(x-x0-rlim))
    else:
        ix1 = -1

    if x[ix0] == x0:
        r = x[ix0:ix1] - x0
    else:
        r = np.insert( x[(ix0+1):ix1] - x0, 0, 0.0)

    f_integral = np.zeros(np.size(r),dtype=np.float32)

    f_interp = interpolate.interp2d(x,y,fxy,kind='cubic')

    Nphi = np.size(x)
    phi_array = np.linspace(0,2*math.pi,num=Nphi)


    for phi in phi_array:
        x_interp = x0 + math.cos(phi) * r
        y_interp = y0 + math.sin(phi) * r

        if 0.0 <= phi < math.pi/2:
            f_integral += np.diagonal(f_interp(x_interp,y_interp))
        elif math.pi/2 <= phi < math.pi:
            f_integral += np.flip(np.diagonal(np.flipud(f_interp(x_interp,y_interp))),axis=0)
        elif math.pi <= phi < 3*math.pi/2:
            f_integral += np.diagonal(np.flipud(np.fliplr(f_interp(x_interp,y_interp))))
        elif 3*math.pi/2 <= phi <= 2*math.pi:
            f_integral += np.flip(np.diagonal(np.fliplr(f_interp(x_interp,y_interp))),axis=0) 

    fr = f_integral/Nphi

    return r, fr

def density_inversion(r_nb, nb, r_psi, psi, zeta_pos, savepath = './plots/equilib', ifplot=True):

    psi_norm = psi - psi[0]

    psi_interp = np.interp(r_nb, r_psi, psi_norm)
    dnb = np.diff(nb)
    dpsi = np.diff(psi_interp)

    x = psi_interp[1:]
    F = -1*np.divide(dnb,dpsi)/(2*math.pi)

    popt, pcov = curve_fit(exp_func, x, F)

    if ifplot:
        mkdirs_if_nexist(savepath)

        fig = plt.figure()
        plt.plot(r_nb, nb)
        ax = plt.gca()        
        ax.set_xlabel(r'$k_p r$', fontsize=14)
        ax.set_ylabel(r'$n_b/n_0$', fontsize=14)   
        saveas_eps_pdf(fig, savepath, 'nb_zeta_%0.2f' % zeta_pos)

        fig = plt.figure()
        plt.plot(r_nb[:-1], np.divide(dnb,np.diff(r_nb)))
        ax = plt.gca()        
        ax.set_xlabel(r'$k_p r$', fontsize=14)
        ax.set_ylabel(r'$d(n_b/n_0)/(k_p dr)$', fontsize=14)   
        saveas_eps_pdf(fig, savepath, 'dnb_dr_zeta_%0.2f' % zeta_pos)

        fig = plt.figure()
        plt.plot(r_nb[:-1], np.divide(dpsi,np.diff(r_nb)))
        ax = plt.gca()        
        ax.set_xlabel(r'$k_p r$', fontsize=14)
        ax.set_ylabel(r'$d(e\psi/mc^2)/(k_p dr)$', fontsize=14)    
        saveas_eps_pdf(fig, savepath, 'dpsi_dr_zeta_%0.2f' % zeta_pos)


        fig = plt.figure()
        plt.plot(psi_norm, r_nb)
        ax = plt.gca()        
        ax.set_ylabel(r'$k_p r$', fontsize=14)
        ax.set_xlabel(r'$e\psi/mc^2$', fontsize=14) 
        saveas_eps_pdf(fig, savepath, 'r_vs_psi_zeta_%0.2f' % zeta_pos)


        fig = plt.figure()
        plt.plot(r_nb, psi_norm)
        ax = plt.gca()        
        ax.set_xlabel(r'$k_p r$', fontsize=14)
        ax.set_ylabel(r'$e\psi/mc^2$', fontsize=14) 
        saveas_eps_pdf(fig, savepath, 'psi_zeta_%0.2f' % zeta_pos)

        fig = plt.figure()
        plt.plot(x, F, '-', label=r'F(x)')
        plt.plot(x, exp_func(x, *popt), '--', label=r'fit: $%0.2f \cdot exp(-%0.2f \cdot x)$' % tuple(popt))
        ax = plt.gca()        
        ax.set_xlabel(r'$x$', fontsize=14)
        ax.set_ylabel(r'$F$', fontsize=14) 
        plt.legend(frameon=False)
        saveas_eps_pdf(fig, savepath, 'F_zeta_%0.2f' % zeta_pos)        

    return x, F


class Equilib:
    def __init__(self, nr, Wr):
        pass


def main():
    parser = parseargs()
    args = parser.parse_args()

    for zeta_pos in args.zeta_pos:
        nbdat = Data_xy_slice(args.nb_path, zeta_pos)
        x_n, y_n, nb_xy = nbdat.read(navg=args.nzetaavg)

        nb_abs_filtered = filters.gaussian_filter(np.abs(nb_xy), args.smooth_sigma, mode='mirror')
        # nb_abs_filtered = np.abs(nb_xy)

        fig = plt.figure()
        cax = plt.pcolor(np.abs(nb_xy))
        ax = plt.gca()        
        ax.set_xlabel(r'$x$', fontsize=14)
        ax.set_ylabel(r'$y$', fontsize=14)
        cbar = fig.colorbar(cax)
        saveas_png(fig, args.savepath, 'nb_%0.2f' % zeta_pos)

        fig = plt.figure()
        cax = plt.pcolor(nb_abs_filtered)
        ax = plt.gca()        
        ax.set_xlabel(r'$x$', fontsize=14)
        ax.set_ylabel(r'$y$', fontsize=14) 
        cbar = fig.colorbar(cax)
        saveas_png(fig, args.savepath, 'nb_filtered_%0.2f' % zeta_pos)

        fig = plt.figure()
        cax = plt.pcolor(np.abs(nb_xy) - nb_abs_filtered)
        ax = plt.gca()        
        ax.set_xlabel(r'$x$', fontsize=14)
        ax.set_ylabel(r'$y$', fontsize=14) 
        cbar = fig.colorbar(cax)
        saveas_png(fig, args.savepath, 'nb_diff_%0.2f' % zeta_pos)

        r_n, nbr = cart_to_r_transform(x_n, y_n, nb_abs_filtered, rlim=args.rlim)

        if (args.mpsi_path != None):
            psidat = Data_xy_slice(args.mpsi_path, zeta_pos)
            x_psi, y_psi, mpsi_xy = psidat.read(navg=args.nzetaavg)

        else:
            wxdat = Data_xy_slice(args.wx_path, zeta_pos)
            x_wx, y_wy, wx_xy  = wxdat.read(navg=args.nzetaavg)
            x_psi, y_psi, mpsi_xy = cart_int(x_wx, y_wy, wx_xy, axis='x')

        r_psi, mpsi = cart_to_r_transform(x_psi, y_psi, mpsi_xy, rlim=args.rlim)

        x, F = density_inversion(r_n, nbr, r_psi, mpsi, zeta_pos = zeta_pos, savepath = args.savepath, ifplot = True)


if __name__ == "__main__":
    main()
