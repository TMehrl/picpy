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
    parser.add_argument(  'zeta_pos', 
                          metavar='ZETA_POS',                          
                          type=float,
                          default=None,                           
                          help = 'Zeta position.')     
    parser.add_argument(  '--wx-path',
                          dest='wx_path',
                          metavar='WX_PATH',
                          type=str,
                          default=None,                           
                          help = 'Path to Wr file.')
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



def cart_int(x,y,fxy,axis='x'):
    if axis == 'x':
        f_int = integrate.cumtrapz(fxy, x, axis=0, initial=0)
    elif axis == 'y':
        f_int = integrate.cumtrapz(fxy, y, axis=1, initial=0)
    else:
        print('ERROR: Only "x" or "y" allowed for axis!')
        sys.exit(1)

    return x, y, f_int

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

    return r, fr

def density_inversion(r_nb, nb, r_psi, psi):
    psi_interp = np.interp(r_nb, r_psi, psi)
    dnb = np.diff(nb)
    dpsi = np.diff(psi_interp)
    dn_per_dpsi = np.divide(dnb,dpsi)

    return psi_interp[:-1], dn_per_dpsi


class Equilib:
    def __init__(self, nr, Wr):
        pass


def main():
    parser = parseargs()
    args = parser.parse_args()

    nbdat = Data_xy_slice(args.nb_path, args.zeta_pos)
    x_n, y_n, nb_xy = nbdat.read()
    r_n, nbr = cart_to_r_transform(x_n, y_n, np.abs(nb_xy),rlim=args.rlim)

    if (args.mpsi_path != None) and (args.wx_path == None):
        psidat = Data_xy_slice(args.mpsi_path, args.zeta_pos)
        x_psi, y_psi, mpsi_xy = psidat.read()

    elif (args.mpsi_path == None) and (args.wx_path != None):
        wxdat = Data_xy_slice(args.wx_path, args.zeta_pos)
        x_wx, y_wy, wx_xy  = wxdat.read()
        x_psi, y_psi, mpsi_xy = cart_int(x_wx, y_wy, wx_xy,axis='x')
    else:
        print('ERROR: Either wx or mpsi path must be specified!')
        sys.exit(1)
    
    r_psi, mpsi = cart_to_r_transform(x_psi, y_psi, mpsi_xy,rlim=args.rlim)

    x, F = density_inversion(r_n, nbr, r_psi, -1*mpsi)

    savepath = './plots/equilib'

    mkdirs_if_nexist(savepath)

    fig = plt.figure()
    plt.plot(r_n, nbr)
    ax = plt.gca()        
    ax.set_xlabel(r'$k_p r$', fontsize=14)
    ax.set_ylabel(r'$n_b/n_0$', fontsize=14)   
    saveas_eps_pdf(fig, savepath, 'nb')

    fig = plt.figure()
    plt.plot(r_psi, mpsi)
    ax = plt.gca()        
    ax.set_xlabel(r'$k_p r$', fontsize=14)
    ax.set_ylabel(r'$e\psi/mc^2$', fontsize=14) 
    saveas_eps_pdf(fig, savepath, 'psi')

    fig = plt.figure()
    plt.plot(mpsi, r_psi)
    ax = plt.gca()        
    ax.set_ylabel(r'$k_p r$', fontsize=14)
    ax.set_xlabel(r'$e\psi/mc^2$', fontsize=14) 
    saveas_eps_pdf(fig, savepath, 'r_vs_psi')

    fig = plt.figure()
    plt.plot(x, F)
    ax = plt.gca()        
    ax.set_xlabel(r'$x$', fontsize=14)
    ax.set_ylabel(r'$F$', fontsize=14) 
    saveas_eps_pdf(fig, savepath, 'F')

if __name__ == "__main__":
    main()
