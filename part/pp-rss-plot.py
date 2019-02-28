#!/usr/bin/env python3
# This script may be executed like this:
# nohup ./raw-slice-series-plotting.py <DATA>/ 1> rss.out 2> rss.err &

import os
import sys
import argparse

import math
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
from matplotlib.ticker import FormatStrFormatter

mypath = os.path.dirname(os.path.abspath( __file__ ))
incpath = os.path.split(mypath)[0] + '/inc'
sys.path.append(incpath)
import pp_defs
from pp_h5dat import SliceMoms
from pp_h5dat import H5Plot
from pp_h5dat import mkdirs_if_nexist
import pp_raw_ana
from pp_plt_tools import saveas_png
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
        prefix = 'rss_name'
        path = './plots'

def ps_parseopts():

    desc="""This is the picpy postprocessing tool."""

    parser = argparse.ArgumentParser(description=desc)
    #parser = OptionParser(usage=usg, description=desc)
    parser.add_argument('path',
                        metavar = 'PATH',
                        help = 'Path to raw files.')    
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        dest="verbose",
                        default=True,
                        help = "Print info (default: %(default)s).")
    parser.add_argument(  "-q", "--quiet",
                        action="store_false",
                        dest="verbose",
                        help = "Don't print info.")
    parser.add_argument("-s", "--save-path",
                        action="store",
                        dest="savepath",
                        metavar="PATH",
                        default=parsedefs.save.path + '/raw-slice-series',
                        help = """Path to which generated files will be saved.
                            (default: %(default)s)""")
    parser.add_argument("-n", "--name-prefix",
                        dest="save_prefix",
                        metavar="NAME",
                        default=parsedefs.save.prefix,
                        help = """Define customized prefix of output filename.""")
    parser.add_argument("-a", "--axis",
                        dest = "axis",
                        metavar="AXIS",
                        choices=['x', 'y', 'r'],
                        default='x',
                        help= """Primary axis for plotting
                              (default: %(default)s).""")
    parser.add_argument("-c", "--code",
                        action='store',
                        dest="piccode",
                        metavar="CODE",
                        choices = [pp_defs.code.hipace, pp_defs.code.osiris,],
                        default = pp_defs.code.hipace,
                        help= "PIC code which was used to generate files (default: %(default)s).")
    parser.add_argument("-d", "--dim",
                        action='store',
                        dest="dimensionality",
                        metavar="DIM",
                        choices=[1, 2, 3,],
                        default=3,
                        help= """Dimensionality of PIC simulation
                              (default: %(default)s).""")
    parser.add_argument("-N", "--number-of-files",
                        action='store',
                        dest="Nfiles",
                        metavar="NFILES",
                        default=0,
                        help= """Number of files to analyze.""")
    parser.add_argument("-o", "--mom-order",
                        type=int,
                        action='store',
                        dest="mom_order",
                        metavar="MOMORDER",
                        choices=[1, 2, 3, 4,],
                        default=None,
                        help='Order of moment evaluation (Default: %(default)s).')    
    parser.add_argument('-t', '--time',
                        help='time for which rms plots are to be generated',
                        action='store',
                        dest="time",
                        nargs=1,
                        type=float,
                        default=None)
    parser.add_argument('--zeta-pos',
                        help='zeta position for which rms plots are to be generated',
                        action='store',
                        dest="zeta_pos",
                        type=float,
                        default=None)    
    parser.add_argument("--use-time",
                        action="store_false",
                        dest="t_is_z",
                        default=True,                          
                        help = "Use time 't' for axes labels (instead of z).")                           
    parser.add_argument('--zeta-range',
                        help='zeta range',
                        action='store',
                        dest="zeta_range",
                        metavar=('ZETA_MIN', 'ZETA_MAX'),
                        nargs=2,
                        type=float,
                        default=None)    
    parser.add_argument("--h5",
                        dest = "h5plot",
                        action="store_true",
                        default=True,
                        help = "Save plot as hdf5 file (Default: %(default)s).") 
    parser.add_argument("--latexon",
                        dest = "latexon",
                        action="store_true",
                        default=False,
                        help = "Use LaTeX font (Default: %(default)s).")
    return parser



def magn_check(x):
    return not (-3.0 < math.log(1e-14 + np.max(abs(x)),10) < 3.0)    


def plot_save_slice_rms(slm, savepath, verbose=True, t_is_z=True):

    x = slm.get_zeta_array()
    y = slm.get_time_array()

    if t_is_z:
        xlabel_str = r'$k_p z$'
    else:
        xlabel_str = r'$\omega_p t$'

    sigma_x = np.sqrt( np.absolute( slm.get(x=2) ) )
    fig_sx = plt.figure()
    cax = plt.pcolormesh( x,
                          y,
                          sigma_x,
                          cmap=cm.Blues,
                          vmin=0, vmax=np.amax(abs(sigma_x)) )
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(xlabel_str, fontsize=14)    
    cbar = fig_sx.colorbar( cax )
    cbar.ax.set_ylabel(r'$k_p \sigma_x$', fontsize=14)
    saveas_png(fig_sx, savepath, 'sigma_x')
    plt.close(fig_sx)

    sigma_px = np.sqrt( np.absolute( slm.get(px=2) ) )
    fig_spx = plt.figure()
    cax = plt.pcolormesh( x,
                          y,
                          sigma_px,
                          cmap=cm.YlGn,
                          vmin=0, vmax=np.amax(sigma_px) )
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(xlabel_str, fontsize=14)    
    cbar = fig_spx.colorbar( cax )
    cbar.ax.set_ylabel(r'$k_p \sigma_{p_x}$', fontsize=14)  
    saveas_png(fig_spx, savepath, 'sigma_px')
    plt.close(fig_spx)

    xpx = slm.get(x=1,px=1)
    fig_xpx = plt.figure()
    cax = plt.pcolormesh( x,
                          y,
                          xpx,
                          cmap=cm.BrBG,
                          vmin=-np.amax(np.abs(xpx)), vmax=np.amax(np.abs(xpx)) )
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(xlabel_str, fontsize=14)    
    cbar = fig_xpx.colorbar( cax )
    cbar.ax.set_ylabel(r'$k_p \left\langle x p_x\right\rangle$', fontsize=14)
    saveas_png(fig_xpx, savepath, 'xpx')
    plt.close(fig_xpx)

    emittance = np.sqrt( np.multiply(slm.get(x=2), slm.get(px=2)) 
                         - np.power(slm.get(x=1,px=1),2) )
    fig_e = plt.figure()
    cax = plt.pcolormesh( x,
                          y,
                          emittance,
                          cmap=cm.Reds,
                          vmin=np.amin(emittance), vmax=np.amax(emittance) )
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(xlabel_str, fontsize=14)    
    cbar = fig_e.colorbar( cax )
    cbar.ax.set_ylabel(r'$k_p \epsilon_x$', fontsize=14)
    saveas_png(fig_e, savepath, 'slice_emittance_x')
    plt.close(fig_e)

def plot_save_proj_rms(slm, savepath, axdir='x', h5plot=True, verbose=True, t_is_z=True):
    
    t = slm.get_time_array()
    tot_charge = np.sum(slm.get_charge(), axis=1)

    if t_is_z:
        xlabel_str = r'$k_p z$'
    else:
        xlabel_str = r'$\omega_p t$'

    if axdir == 'x':
        xsq = slm.get_proj(x=2)
        psq = slm.get_proj(px=2)
        xp = slm.get_proj(x=1,px=1)
        
        emittance_all_slices = np.sqrt( np.multiply(slm.get(x=2), slm.get(px=2)) - np.power(slm.get(x=1,px=1),2) )
        
        sigma_xyr_lab = r'$k_p\sigma_x$'
        sigma_xyr_savename = 'sigma_x_proj'
        
        sigma_pxyr_lab = r'$\sigma_{p_x}/m_e c$'
        sigma_pxyr_savename = 'sigma_px_proj'
        
        corr_xy_lab = r'$ k_p \left \langle x\,p_x \right\rangle/m_e c$'
        corr_xy_savename = 'xpx_proj'
        
        epsilon_xy_proj_lab = r'$k_p \epsilon_x$'
        epsilon_xy_proj_savename = 'emittance_x_proj'   

        epsilon_xy_sl_lab = r'$k_p \epsilon_{x,\mathrm{sliced}}$'
        epsilon_xy_sl_savename = 'emittance_x_sliced' 


    elif axdir == 'y':
        xsq = slm.get_proj(y=2)
        psq = slm.get_proj(py=2)
        xp = slm.get_proj(y=1,py=1)
        
        emittance_all_slices = np.sqrt( np.multiply(slm.get(y=2), slm.get(py=2)) - np.power(slm.get(y=1,py=1),2) )
        
        sigma_xyr_lab = r'$k_p\sigma_y$'
        sigma_xyr_savename = 'sigma_y_proj'                
        
        sigma_pxyr_lab = r'$\sigma_{p_y}/m_e c$'
        sigma_pxyr_savename = 'sigma_py_proj'
        
        corr_xy_lab = r'$ k_p \left \langle y\,p_y \right\rangle/m_e c$'
        corr_xy_savename = 'ypy_proj'

        epsilon_xy_proj_lab = r'$k_p \epsilon_y$'
        epsilon_xy_proj_savename = 'emittance_y_proj'

        epsilon_xy_sl_lab = r'$k_p \epsilon_{y,\mathrm{sliced}}$'
        epsilon_xy_sl_savename = 'emittance_y_sliced'         

    elif axdir == 'r':
        xavg = slm.project(slm.get(x=1) + slm.get(y=1))
        pavg = slm.project(slm.get(px=1) + slm.get(py=1))

        xsq_noncentral = slm.project(slm.get(x=2) + slm.get(y=2) \
            + np.power(slm.get(x=1),2) + np.power(slm.get(y=1),2))
        psq_noncentral = slm.project(slm.get(px=2) \
            + slm.get(py=2) + np.power(slm.get(px=1),2) + np.power(slm.get(py=1),2))
        xp_noncentral = slm.project(slm.get(x=1,px=1) + slm.get(y=1,py=1) \
            + np.multiply(slm.get(x=1),slm.get(px=1)) + np.multiply(slm.get(y=1),slm.get(py=1)))

        xsq = xsq_noncentral - np.pow(xavg,2)
        psq = psq_noncentral - np.pow(pavg,2)
        xp = xp_noncentral = np.multiply(xavg,pavg)

        emittance_all_slices = np.sqrt( np.multiply(slm.get(x=2), slm.get(px=2)) - np.power(slm.get(x=1,px=1),2) + \
                np.multiply(slm.get(y=2), slm.get(py=2)) - np.power(slm.get(y=1,py=1),2) )
        
        sigma_xyr_lab = r'$k_p\sigma_r$'
        sigma_xyr_savename = 'sigma_r_proj'                
        
        sigma_pxyr_lab = r'$(\sigma_{p_x}^2 + \sigma_{p_y}^2)^{1/2}/m_e c$'
        sigma_pxyr_savename = 'sigma_pr_proj'
        
        corr_xy_lab = r'$ k_p (\left \langle x\,p_x \right\rangle + \left \langle y\,p_y \right\rangle)/m_e c$'
        corr_xy_savename = 'xpx_ypy_proj'

        epsilon_xy_proj_lab = r'$k_p \epsilon_r$'
        epsilon_xy_proj_savename = 'emittance_r_proj'

        epsilon_xy_sl_lab = r'$k_p \epsilon_{r,\mathrm{sliced}}$'
        epsilon_xy_sl_savename = 'emittance_r_sliced' 


    emittance_sliced = slm.project(emittance_all_slices)

    fig_sx = plt.figure()    
    plt.plot(t, np.sqrt(xsq))
    ax = plt.gca()
    ax.set_xlabel(xlabel_str, fontsize=14) 
    ax.set_ylabel(sigma_xyr_lab, fontsize=14)
    if magn_check(np.sqrt(xsq)):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15) 
    
    saveas_eps_pdf(fig_sx, savepath, sigma_xyr_savename, h5plot=h5plot)                  
    plt.close(fig_sx)

    fig_sp = plt.figure()    
    plt.plot(t, np.sqrt(psq))
    ax = plt.gca()
    ax.set_xlabel(xlabel_str, fontsize=14) 
    ax.set_ylabel(sigma_pxyr_lab, fontsize=14)
    if magn_check(np.sqrt(psq)):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15)  
    
    saveas_eps_pdf(fig_sp, savepath, sigma_pxyr_savename, h5plot=h5plot)                  
    plt.close(fig_sp)


    fig_xp = plt.figure()    
    plt.plot(t, xp)
    ax = plt.gca()
    ax.set_xlabel(xlabel_str, fontsize=14) 
    ax.set_ylabel(corr_xy_lab, fontsize=14)
    if magn_check(xp):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15)
    
    saveas_eps_pdf(fig_xp, savepath, corr_xy_savename, h5plot=h5plot)   
    plt.close(fig_xp)


    emittance_proj = np.sqrt(np.multiply(xsq,psq)-np.power(xp,2))
    fig_e = plt.figure()    
    plt.plot(t, emittance_proj, label='projected emittance') 
    ax = plt.gca()
    ax.set_xlabel(xlabel_str, fontsize=14) 
    ax.set_ylabel(epsilon_xy_proj_lab, fontsize=14)
    if magn_check(emittance_proj):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15)  
    
    saveas_eps_pdf(fig_e, savepath, epsilon_xy_proj_savename, h5plot=h5plot)   
    plt.close(fig_e)

    fig_esl = plt.figure()    
    plt.plot(t, emittance_sliced, label='sliced emittance') 
    ax = plt.gca()
    ax.set_xlabel(xlabel_str, fontsize=14) 
    ax.set_ylabel(epsilon_xy_sl_lab, fontsize=14)
    if magn_check(emittance_sliced):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15)  
    
    saveas_eps_pdf(fig_esl, savepath, epsilon_xy_sl_savename, h5plot=h5plot)   
    plt.close(fig_esl)

    fig_e_slpr = plt.figure()    
    epp = plt.plot(t, emittance_proj, label='projected')
    esp = plt.plot(t, emittance_sliced, color = epp[0].get_color(), linestyle='--', label='sliced')    
    ax = plt.gca()
    ax.set_xlabel(xlabel_str, fontsize=14) 
    ax.set_ylabel(epsilon_xy_proj_lab, fontsize=14)
    ax.legend()
    if magn_check(emittance_proj):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15)  
    
    saveas_eps_pdf(fig_e_slpr, savepath, epsilon_xy_proj_savename + '_sl', h5plot=h5plot)   
    plt.close(fig_e_slpr)


def plot_save_slice_rms_lines(slm, savepath, time = None, axdir='x', h5plot=True):

    if time == None:
        tidx = [0,-1];
    else:
        tidx = [(np.abs(slm.get_time_array() - time)).argmin()]

    for i in tidx:
        if axdir == 'x':
            sigma_xy = np.sqrt(slm.get(x=2)[i,:])
            sigma_xy_lab = r'$k_p \sigma_{x}$'
            sigma_xy_savename = 'sigma_x'
            sigma_pxy = np.sqrt(slm.get(px=2)[i,:])
            sigma_pxy_lab = r'$\sigma_{p_x}/mc$'
            sigma_pxy_savename = 'sigma_px'
            # Also define labels and savenames here!    
        elif axdir == 'y':
            sigma_xy = np.sqrt(slm.get(y=2)[i,:])
            sigma_xy_lab = r'$k_p \sigma_{y}$'
            sigma_xy_savename = 'sigma_y'
            sigma_pxy = np.sqrt(slm.get(py=2)[i,:])
            sigma_pxy_lab = r'$\sigma_{p_y}/mc$'
            sigma_pxy_savename = 'sigma_py'
        elif axdir == 'r':
            sigma_xy = np.sqrt(slm.get(x=2)[i,:] + slm.get(y=2)[i,:])
            sigma_xy_lab = r'$k_p \sigma_{r}$'
            sigma_xy_savename = 'sigma_r'
            sigma_pxy = np.sqrt(slm.get(px=2)[i,:] + slm.get(py=2)[i,:])
            sigma_pxy_lab = r'$(\sigma_{p_x}^2+\sigma_{p_y}^2)^{1/2}/mc$'
            sigma_pxy_savename = 'sigma_pr'            

        fig_sigma_xy = plt.figure()
        plt.plot( slm.get_zeta_array(),
                  sigma_xy)
        ax = plt.gca()
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
        ax.set_ylabel(sigma_xy_lab, fontsize=14)     
        saveas_eps_pdf(fig_sigma_xy, savepath, ('%s_time_%0.1f' % (sigma_xy_savename, slm.get_time_array()[i])), h5plot=h5plot)
        plt.close(fig_sigma_xy)


        fig_sigma_pxy = plt.figure()
        plt.plot( slm.get_zeta_array(),
                  sigma_pxy)
        ax = plt.gca()
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
        ax.set_ylabel(sigma_pxy_lab, fontsize=14)     
        saveas_eps_pdf(fig_sigma_pxy, savepath, ('%s_time_%0.1f' % (sigma_pxy_savename, slm.get_time_array()[i])), h5plot=h5plot)
        plt.close(fig_sigma_pxy)

def plot_save_slice_exkurtosis_lines(slm, savepath, time = None, axdir='x', h5plot=True):

    if time == None:
        tidx = [0,-1];
    else:
        tidx = [(np.abs(slm.get_time_array() - time)).argmin()]

    for i in tidx:
        if axdir == 'y':
            var_xy = slm.get(y=2)[i,:]
            var_xy_nonzero_idx = np.nonzero(var_xy)[0]           
            exkurtosis_xy = np.divide(  slm.get(y=4)[i,var_xy_nonzero_idx],\
                                        np.power(var_xy[var_xy_nonzero_idx],2)) - 3
            exkurtosis_xy_lab = r'$\left \langle y^4 \right\rangle -3$'
            exkurtosis_xy_savename = 'exkurtosis_y'
            
            var_pxy = slm.get(py=2)[i,:]
            var_pxy_nonzero_idx = np.nonzero(var_pxy)[0]
            exkurtosis_pxy = np.divide( slm.get(py=4)[i,var_pxy_nonzero_idx],\
                                        np.power(var_pxy[var_pxy_nonzero_idx],2)) - 3
            exkurtosis_pxy_lab = r'$\left \langle p_y^4/\sigma_{py}^4 \right\rangle -3$'
            exkurtosis_pxy_savename = 'exkurtosis_py'

        else: 
            # xdir == 'x' or xdir == 'r'
            var_xy = slm.get(x=2)[i,:]
            var_xy_nonzero_idx = np.nonzero(var_xy)[0]
            exkurtosis_xy = np.divide(  slm.get(x=4)[i,var_xy_nonzero_idx],\
                                        np.power(var_xy[var_xy_nonzero_idx],2)) - 3
            exkurtosis_xy_lab = r'$\left \langle x^4/\sigma_x^4 \right\rangle - 3$'
            exkurtosis_xy_savename = 'exkurtosis_x'
            
            var_pxy = slm.get(px=2)[i,:]
            var_pxy_nonzero_idx = np.nonzero(var_pxy)[0]
            exkurtosis_pxy = np.divide( slm.get(px=4)[i,var_pxy_nonzero_idx],\
                                        np.power(var_pxy[var_pxy_nonzero_idx],2)) - 3
            exkurtosis_pxy_lab = r'$\left \langle p_x^4/\sigma_{px}^4 \right\rangle -3$'
            exkurtosis_pxy_savename = 'exkurtosis_px'
            # Also define labels and savenames here!   

        fig_exkurtosis_xy = plt.figure()
        plt.plot( slm.get_zeta_array()[var_xy_nonzero_idx],
                  exkurtosis_xy)
        ax = plt.gca()
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
        ax.set_ylabel(exkurtosis_xy_lab, fontsize=14)     
        saveas_eps_pdf(fig_exkurtosis_xy, savepath, ('%s_time_%0.1f' % (exkurtosis_xy_savename, slm.get_time_array()[i])), h5plot=h5plot)
        plt.close(fig_exkurtosis_xy)


        fig_exkurtosis_pxy = plt.figure()
        plt.plot( slm.get_zeta_array()[var_pxy_nonzero_idx],
                  exkurtosis_pxy)
        ax = plt.gca()
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
        ax.set_ylabel(exkurtosis_pxy_lab, fontsize=14)     
        saveas_eps_pdf(fig_exkurtosis_pxy, savepath, ('%s_time_%0.1f' % (exkurtosis_pxy_savename, slm.get_time_array()[i])), h5plot=h5plot)
        plt.close(fig_exkurtosis_pxy)


def plot_save_slice_quad_corr_lines(slm, savepath, time = None, axdir='x', h5plot=True):

    if time == None:
        tidx = [0,-1];
    else:
        tidx = [(np.abs(slm.get_time_array() - time)).argmin()]

    for i in tidx:
        if axdir == 'y':
            var_xy = slm.get(y=2)[i,:]
            var_pxy = slm.get(py=2)[i,:]
            nonzero_idx = np.nonzero(np.multiply(var_xy, var_pxy))[0]
            quad_corr_xy = np.divide(  slm.get(y=2,py=2)[i,nonzero_idx],\
                                        np.multiply(var_xy[nonzero_idx],var_pxy[nonzero_idx]) ) - 1.0
            quad_corr_xy_lab = r'$\left \langle y^2  p_y^2\right\rangle/(\sigma_y^2\sigma_{p_y}^2) - 1$'
            quad_corr_xy_savename = 'quad_corr_y'
  
        else: 
            # axdir == 'x' or axdir == 'r'
            var_xy = slm.get(x=2)[i,:]
            var_pxy = slm.get(px=2)[i,:]
            nonzero_idx = np.nonzero(np.multiply(var_xy, var_pxy))[0]
            quad_corr_xy = np.divide(  slm.get(x=2,px=2)[i,nonzero_idx],\
                                        np.multiply(var_xy[nonzero_idx],var_pxy[nonzero_idx]) ) - 1.0
            quad_corr_xy_lab = r'$\left \langle x^2  p_x^2\right\rangle/(\sigma_x^2\sigma_{p_x}^2) - 1$'
            quad_corr_xy_savename = 'quad_corr_x'


        fig_quad_corr_xy = plt.figure()
        plt.plot( slm.get_zeta_array()[nonzero_idx],
                  quad_corr_xy)
        ax = plt.gca()
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
        ax.set_ylabel(quad_corr_xy_lab, fontsize=14)     
        saveas_eps_pdf(fig_quad_corr_xy, savepath, ('%s_time_%0.1f' % (quad_corr_xy_savename, slm.get_time_array()[i])), h5plot=h5plot)
        plt.close(fig_quad_corr_xy)


def plot_save_slice_centroids(slm, savepath, h5plot=True, time=None, zeta_pos=None, t_is_z=True):

    if t_is_z:
        xlabel_str = r'$k_p z$'
    else:
        xlabel_str = r'$\omega_p t$'


    Xb0_for_norm = np.ones(slm.get(x=1)[0,:].shape)
    Yb0_for_norm = np.ones(slm.get(y=1)[0,:].shape)
    for i in range(0,len(slm.get_zeta_array())):
        if (slm.get(x=1)[0,i] != 0.0):
            Xb0_for_norm[i] = slm.get(x=1)[0,i]
            Yb0_for_norm[i] = slm.get(y=1)[0,i]

    if time == None:
        time = slm.get_time_array()[-1]
        tidx = [-1];
    else:
        tidx = [(np.abs(slm.get_time_array() - time)).argmin()]

    if len(tidx) == 1:
        time = [time]

    for i in range(0,len(tidx)):
        Xb_lout_savename = 'Xb_time_%0.2f' % time[i]
        Yb_lout_savename = 'Yb_time_%0.2f' % time[i] 
        Xb_lout_ylab = r'$k_p X_{b}(t=%0.1f)$' % time[i]
        Yb_lout_ylab = r'$k_p Y_{b}(t=%0.1f)$' % time[i]

        figXb_lout = plt.figure()
        plt.plot(slm.get_zeta_array(), slm.get(x=1)[tidx[i],:])
        ax = plt.gca()
        ymin, ymax = ax.get_ylim()
        if ymin > 0 and ymax > 0:
            plt.ylim(0, ymax*1.2)
        elif ymin < 0 and ymax < 0:
            plt.ylim(ymin*1.2, 0)        
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
        ax.set_ylabel(Xb_lout_ylab, fontsize=14)
        saveas_eps_pdf(figXb_lout, savepath, Xb_lout_savename, h5plot=h5plot)
        plt.close(figXb_lout)

        figYb_lout = plt.figure()
        plt.plot(slm.get_zeta_array(), slm.get(y=1)[tidx[i],:])
        ax = plt.gca()
        ymin, ymax = ax.get_ylim()
        if ymin > 0 and ymax > 0:
            plt.ylim(0, ymax*1.2)
        elif ymin < 0 and ymax < 0:
            plt.ylim(ymin*1.2, 0)    
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
        ax.set_ylabel(Yb_lout_ylab, fontsize=14)  
        saveas_eps_pdf(figYb_lout, savepath, Yb_lout_savename, h5plot=h5plot)
        plt.close(figYb_lout)

    Xb_norm = np.zeros( slm.get(x=1).shape )
    Yb_norm = np.zeros( slm.get(y=1).shape )
#    zeta_hseed = np.min(slm.get_zeta_array())
#    idx_hseed = (np.abs(slm.get_zeta_array()-zeta_hseed)).argmin()

    for i in range(0,len(slm.get_zeta_array())):
#        if (slm.get_zeta_array()[i] <= zeta_hseed):
        Xb_norm[:,i] = np.absolute( slm.get(x=1)[:,i]/Xb0_for_norm[i] )
        Yb_norm[:,i] = np.absolute( slm.get(y=1)[:,i]/Yb0_for_norm[i] )

    figXb = plt.figure()
    cax = plt.pcolormesh(   slm.get_zeta_array(),
                            slm.get_time_array(),
                            slm.get(x=1),
                            cmap=cm.PuOr,
                            vmin=-np.amax(abs(slm.get(x=1))), 
                            vmax=np.amax(abs(slm.get(x=1))))
    cbar = figXb.colorbar(cax)   
    cbar.ax.set_ylabel('$k_p X_b$', fontsize=14)
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(xlabel_str, fontsize=14) 
    saveas_png(figXb, savepath, 'Xb')
    plt.close(figXb)

    figYb = plt.figure()
    cax = plt.pcolormesh(   slm.get_zeta_array(),
                            slm.get_time_array(),
                            slm.get(y=1),
                            cmap=cm.PuOr,
                            vmin=-np.amax(abs(slm.get(y=1))), 
                            vmax=np.amax(abs(slm.get(y=1))))
    cbar = figYb.colorbar(cax)
    cbar.ax.set_ylabel('$k_p Y_b$', fontsize=14)
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(xlabel_str, fontsize=14)     
    saveas_png(figYb, savepath, 'Yb')
    plt.close(figYb)

    figXbnorm = plt.figure()
    cax = plt.pcolormesh(   slm.get_zeta_array(),
                            slm.get_time_array(),
                            Xb_norm,
                            cmap=cm.Blues,
                            vmin=0, 
                            vmax=np.amax(Xb_norm))
    cbar = figXbnorm.colorbar(cax)
    cbar.ax.set_ylabel('$|X_b/X_{b,0}|$', fontsize=14)
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(xlabel_str, fontsize=14)  
    saveas_png(figXbnorm, savepath, 'Xb_rel')
    plt.close(figXbnorm)

    figYbnorm = plt.figure()
    cax = plt.pcolormesh(   slm.get_zeta_array(),
                            slm.get_time_array(),
                            Yb_norm,
                            cmap=cm.Blues,
                            vmin=0, 
                            vmax=np.amax(Yb_norm))
    cbar = figYbnorm.colorbar(cax)
    cbar.ax.set_ylabel('$|Y_b/Y_{b,0}|$', fontsize=14)
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(xlabel_str, fontsize=14)
    saveas_png(figYbnorm, savepath, 'Yb_rel')
    plt.close(figYbnorm)

    figXb0 = plt.figure()
    plt.plot(slm.get_zeta_array(), slm.get(x=1)[0,:])
    ax = plt.gca()
    ymin, ymax = ax.get_ylim()
    if ymin > 0 and ymax > 0:
        plt.ylim(0, ymax*1.2)
    elif ymin < 0 and ymax < 0:
        plt.ylim(ymin*1.2, 0)        
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$k_p X_{b,0}$', fontsize=14)
    saveas_eps_pdf(figXb0, savepath, 'Xb0', h5plot=h5plot)
    plt.close(figXb0)

    figYb0 = plt.figure()
    plt.plot(slm.get_zeta_array(), slm.get(y=1)[0,:])
    ax = plt.gca()
    ymin, ymax = ax.get_ylim()
    if ymin > 0 and ymax > 0:
        plt.ylim(0, ymax*1.2)
    elif ymin < 0 and ymax < 0:
        plt.ylim(ymin*1.2, 0)    
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$k_p Y_{b,0}$', fontsize=14)     
    saveas_eps_pdf(figYb0, savepath, 'Yb0', h5plot=h5plot)
    plt.close(figYb0)

    if zeta_pos != None:
        zetaidx = [(np.abs(slm.get_zeta_array() - zeta_pos)).argmin()]
        Xb_savename = 'Xb_zeta_%0.2f' % zeta_pos
        Yb_savename = 'Yb_zeta_%0.2f' % zeta_pos 
    else:
        zetaidx = 0
        Xb_savename = 'Xb_tail'
        Yb_savename = 'Yb_tail'

    figXbtail = plt.figure()
    plt.plot(slm.get_time_array(), slm.get(x=1)[:,zetaidx])
    ax = plt.gca()
    ax.set_xlabel(xlabel_str, fontsize=14)
    ax.set_ylabel(r'$k_p X_{b,\mathrm{tail}}$', fontsize=14)
    if magn_check(slm.get(x=1)[:,zetaidx]):    
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15)
    saveas_eps_pdf(figXbtail, savepath, Xb_savename, h5plot=h5plot)    
    plt.close(figXbtail)

    figYbtail = plt.figure()
    plt.plot(slm.get_time_array(), slm.get(y=1)[:,zetaidx])
    ax = plt.gca()
    ax.set_xlabel(xlabel_str, fontsize=14)
    ax.set_ylabel(r'$k_p Y_{b,\mathrm{tail}}$', fontsize=14)
    if magn_check(slm.get(y=1)[:,zetaidx]):    
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15)    
    saveas_eps_pdf(figYbtail, savepath, Yb_savename, h5plot=h5plot)   
    plt.close(figYbtail)  



def plot_save_slice_ene(slm, savepath, time = None, h5plot=True, t_is_z=True):

    if t_is_z:
        xlabel_str = r'$k_p z$'
    else:
        xlabel_str = r'$\omega_p t$'
    
    gamma = np.sqrt( 1 + np.power(slm.get(pz=1),2) 
                       + np.power(slm.get(px=1),2)
                       + np.power(slm.get(py=1),2) )

    if time == None:
        tidx = [0,-1];
    else:
        tidx = [(np.abs(slm.get_time_array() - time)).argmin()]

    for i in tidx:
        fig = plt.figure()
        plt.plot(slm.get_zeta_array(), gamma[i,:])
        ax = plt.gca()
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
        ax.set_ylabel(r'$\gamma$', fontsize=14)
        saveas_eps_pdf(fig, savepath, ('gamma_time_%0.1f' % (slm.get_time_array()[i])), h5plot=h5plot)
        plt.close(fig)    



    figG = plt.figure()
    cax = plt.pcolormesh( slm.get_zeta_array(),
                          slm.get_time_array(),
                          gamma,
                          cmap=cm.GnBu,
                          vmin=np.amin(gamma), vmax=np.amax(gamma))
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(xlabel_str, fontsize=14)    
    cbar = figG.colorbar( cax )
    cbar.ax.set_ylabel(r'$\gamma$', fontsize=14)
    saveas_png(figG, savepath, 'gamma')
    plt.close(figG)    


def plot_save_ene_ene_spread(slm, savepath, h5plot=True, t_is_z=True):

    if t_is_z:
        xlabel_str = r'$k_p z$'
    else:
        xlabel_str = r'$\omega_p t$'

    t = slm.get_time_array()
    tot_charge = np.sum(slm.get_charge(), axis=1)
    
    avg_gamma = np.divide(np.sum( np.multiply( slm.get(pz=1), slm.get_charge()), axis=1),tot_charge)
    sigma_gamma = np.divide(np.sqrt( np.abs( np.divide(np.sum( np.multiply( slm.get(pz=2) \
        + np.power(slm.get(pz=1),2), slm.get_charge()), axis=1),tot_charge) - np.power(avg_gamma,2))),avg_gamma)

    fig_sg = plt.figure()    
    plt.plot(t, sigma_gamma)
    ax = plt.gca()
    if magn_check(sigma_gamma):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15)     
    ax.set_xlabel(xlabel_str, fontsize=14) 
    ax.set_ylabel(r'$\sigma_\gamma/\gamma$', fontsize=14)  
    saveas_eps_pdf(fig_sg, savepath, 'sigma_gamma', h5plot=h5plot)                  
    plt.close(fig_sg)


    fig_g = plt.figure()    
    plt.plot(t, avg_gamma)
    ax = plt.gca()
    ax.set_xlabel(xlabel_str, fontsize=14) 
    ax.set_ylabel(r'$\gamma$', fontsize=14)    
    saveas_eps_pdf(fig_g, savepath, 'avg_gamma', h5plot=h5plot)                  
    plt.close(fig_g)


def plot_curr_profile(slm, savepath, time = None, h5plot=True, t_is_z=True):
    dzeta = abs(slm.get_zeta_array()[1] - slm.get_zeta_array()[0]);
    curr = slm.get_charge() / dzeta

    if t_is_z:
        xlabel_str = r'$k_p z$'
    else:
        xlabel_str = r'$\omega_p t$'

    # print('Q = %0.3e' % (np.sum(curr) * dzeta ))

    # I_A = 4 * pi * epsilon_0 * m * c^3 / e
    # [curr] = epsilon_0 * m * c^3 / e
    # curr * 4 * pi = I_b/I_A
    Ib_per_IA = curr/(4 * math.pi)

    if time == None:
        tidx = [0,-1];
    else:
        tidx = [(np.abs(slm.get_time_array() - time)).argmin()]

    for i in tidx:
        fig = plt.figure()
        plt.plot(slm.get_zeta_array(), Ib_per_IA[i,:])
        ax = plt.gca()
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
        ax.set_ylabel(r'$I_{b}/I_A$', fontsize=14)
        saveas_eps_pdf(fig, savepath, ('Ib_time_%0.1f' % (slm.get_time_array()[i])), h5plot=h5plot)
        plt.close(fig)    


    figIb = plt.figure()
    cax = plt.pcolormesh(   slm.get_zeta_array(),
                            slm.get_time_array(),
                            Ib_per_IA,
                            cmap=cm.Greys,
                            vmin=0, 
                            vmax=np.amax(abs(Ib_per_IA)))
    cbar = figIb.colorbar(cax)   
    cbar.ax.set_ylabel('$I_b/I_A$', fontsize=14)
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(xlabel_str, fontsize=14) 
    saveas_png(figIb, savepath, 'Ib')
    plt.close(figIb)


def plot_xterms(slm, savepath, time = None, h5plot=True, t_is_z=True):
    dzeta = abs(slm.get_zeta_array()[1] - slm.get_zeta_array()[0]);
    curr = slm.get_charge() / dzeta

    if t_is_z:
        xlabel_str = r'$k_p z$'
    else:
        xlabel_str = r'$\omega_p t$'

    if time == None:
        tidx = [0,-1];
    else:
        tidx = [(np.abs(slm.get_time_array() - time)).argmin()]

    for i in tidx:
        fig = plt.figure()
        plt.plot(slm.get_zeta_array(), slm.get(x=1,y=1)[i,:])
        ax = plt.gca()
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
        ax.set_ylabel(r'$<xy>$', fontsize=14)
        saveas_eps_pdf(fig, savepath, ('xy_%0.1f' % (slm.get_time_array()[i])), h5plot=h5plot)
        plt.close(fig)    



def main():

    parser = ps_parseopts()

    args = parser.parse_args()

    file = args.path

    slm = SliceMoms()
    slm.read(file)

    if slm.get_h5order(file) > 0:
        if args.mom_order != None:
            mom_order = args.mom_order
        else:
            mom_order = slm.get_h5order(file)
    else:
        print('Error:\tMoment order in file <= 0!')
        sys.exit(1)                 

    if args.zeta_range != None:
        slm.truncate_zeta_region(args.zeta_range[0], args.zeta_range[1])

    mkdirs_if_nexist(args.savepath)

    if args.latexon:
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif') 

    plot_curr_profile(slm, args.savepath, time = args.time, h5plot=args.h5plot)

    plot_save_slice_ene(slm, args.savepath, h5plot=args.h5plot)

    plot_save_slice_centroids(slm, args.savepath, h5plot=args.h5plot, zeta_pos=args.zeta_pos)
    
    if mom_order > 1:
        plot_save_slice_rms(slm, args.savepath)
        plot_save_proj_rms(slm, args.savepath, axdir=args.axis, h5plot=args.h5plot)        
        plot_save_slice_rms_lines(slm, args.savepath, time = args.time, axdir=args.axis, h5plot=args.h5plot)
        plot_save_ene_ene_spread(slm, args.savepath, t_is_z=args.t_is_z)

    if slm.get_if_xterms():
        plot_xterms(slm, args.savepath, time = args.time, h5plot=args.h5plot)

    if mom_order > 3:
        plot_save_slice_exkurtosis_lines(slm, args.savepath, time = args.time, axdir=args.axis, h5plot=args.h5plot)
        plot_save_slice_quad_corr_lines(slm, args.savepath, time = args.time, axdir=args.axis, h5plot=args.h5plot)

if __name__ == "__main__":
    main()
