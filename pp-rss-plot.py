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
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib.colors import LogNorm
from matplotlib.ticker import FormatStrFormatter
import pp_defs
import h5py
from pp_h5dat import SliceMoms
from pp_h5dat import H5Plot
from pp_h5dat import mkdirs_if_nexist
import pp_raw_ana


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
    parser.add_argument(  'path',
                          metavar = 'PATH',
                          help = 'Path to raw files.')    
    parser.add_argument(  "-v", "--verbose",
                        action="store_true",
                        dest="verbose",
                        default=True,
                        help = "Print info (default: %(default)s).")
    parser.add_argument(  "-q", "--quiet",
                        action="store_false",
                        dest="verbose",
                        help = "Don't print info.")
    parser.add_argument(  "-s", "--save-path",
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
    parser.add_argument(  '-t', '--time',
                          help='time for which rms plots are to be generated',
                          action='store',
                          dest="time",
                          nargs=1,
                          type=float,
                          default=None)  
    parser.add_argument(  '--zeta-range',
                          help='zeta range',
                          action='store',
                          dest="zeta_range",
                          metavar=('ZETA_MIN', 'ZETA_MAX'),
                          nargs=2,
                          type=float,
                          default=None)    
    parser.add_argument(  "--h5",
                          dest = "h5plot",
                          action="store_true",
                          default=True,
                          help = "Save plot as hdf5 file (Default: %(default)s).") 

    return parser


def saveas_eps_pdf(fig, savepath, savename, h5plot=True, verbose=True, fformat='pdf'):

    fformat_list = ['eps','pdf']

    if fformat not in fformat_list:
        print("Error: fformat must be one of: " + fformat_list)

    fig.savefig( savepath + '/' + savename + '.' + fformat,
                  format=fformat)
      
    if verbose: 
        sys.stdout.write('Saved "%s.%s" at: %s/\n' % (savename, fformat, savepath))
        sys.stdout.flush()

    if h5plot: 
        h5lp = H5Plot()
        h5lp.inherit_matplotlib_line_plots(plt.gca())    
        h5lp.write(savepath + '/' + savename + '.h5')

        if verbose: 
            sys.stdout.write('Saved "%s.h5" at: %s/\n' % (savename, savepath))
            sys.stdout.flush()  


def saveas_png(fig, savepath, savename, verbose=True):

    fformat = 'png'
    fig.savefig( savepath + '/' + savename + '.' + fformat,
                  format=fformat,
                  dpi=600)    
      
    if verbose: 
        sys.stdout.write('Saved "%s.%s" at: %s/\n' % (savename, fformat, savepath))
        sys.stdout.flush()  


def magn_check(x):
    return not (-3.0 < math.log(1e-14 + np.max(abs(x)),10) < 3.0)    


def plot_save_slice_rms(slm, savepath, verbose=True):

    x = slm.zeta_array
    y = slm.time_array

    sigma_x = np.sqrt( np.absolute( slm.avgx2sq ) )
    fig_sx = plt.figure()
    cax = plt.pcolormesh( x,
                          y,
                          sigma_x,
                          cmap=cm.Blues,
                          vmin=0, vmax=np.amax(abs(sigma_x)) )
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$\omega_p t$', fontsize=14)    
    cbar = fig_sx.colorbar( cax )
    cbar.ax.set_ylabel(r'$k_p \sigma_x$', fontsize=14)
    saveas_png(fig_sx, savepath, 'sigma_x')


    sigma_px = np.sqrt( np.absolute( slm.avgp2sq ) )
    fig_spx = plt.figure()
    cax = plt.pcolormesh( x,
                          y,
                          sigma_px,
                          cmap=cm.YlGn,
                          vmin=0, vmax=np.amax(sigma_px) )
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$\omega_p t$', fontsize=14)    
    cbar = fig_spx.colorbar( cax )
    cbar.ax.set_ylabel(r'$k_p \sigma_{p_x}$', fontsize=14)  
    saveas_png(fig_spx, savepath, 'sigma_px')

    emittance = np.sqrt( np.multiply(slm.avgx2sq, slm.avgp2sq) 
                         - np.power(slm.avgx2p2,2) )
    fig_e = plt.figure()
    cax = plt.pcolormesh( x,
                          y,
                          emittance,
                          cmap=cm.Reds,
                          vmin=np.amin(emittance), vmax=np.amax(emittance) )
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$\omega_p t$', fontsize=14)    
    cbar = fig_e.colorbar( cax )
    cbar.ax.set_ylabel(r'$k_p \epsilon_x$', fontsize=14)
    saveas_png(fig_e, savepath, 'slice_emittance_x')


def plot_save_proj_rms(slm, savepath, axdir=2, h5plot=True, verbose=True):
    
    t = slm.time_array
    tot_charge = np.sum(slm.charge, axis=1)

    if axdir == 2:
        xsq = np.divide(np.sum(np.multiply(slm.avgx2sq, slm.charge), axis=1),tot_charge)
        psq = np.divide(np.sum(np.multiply(slm.avgp2sq, slm.charge), axis=1),tot_charge)
        xp = np.divide(np.sum(np.multiply(slm.avgx2p2, slm.charge), axis=1),tot_charge)
        emittance_all_slices = np.sqrt( np.multiply(slm.avgx2sq, slm.avgp2sq) - np.power(slm.avgx2p2,2) )
        emittance_sliced = np.divide(np.sum(np.multiply(emittance_all_slices, slm.charge), axis=1),tot_charge)
        # Also define labels and savenames here!    
    elif axdir == 3:
        xsq = np.divide(np.sum(np.multiply(slm.avgx3sq, slm.charge), axis=1),tot_charge)
        psq = np.divide(np.sum(np.multiply(slm.avgp3sq, slm.charge), axis=1),tot_charge)
        xp = np.divide(np.sum(np.multiply(slm.avgx3p3, slm.charge), axis=1),tot_charge)
        emittance_all_slices = np.sqrt( np.multiply(slm.avgx3sq, slm.avgp3sq) - np.power(slm.avgx3p3,2) )
        emittance_sliced = np.divide(np.sum(np.multiply(emittance_all_slices, slm.charge), axis=1),tot_charge)
        # Also define labels and savenames here!           

    fig_sx = plt.figure()    
    plt.plot(t, np.sqrt(xsq))
    ax = plt.gca()
    ax.set_xlabel(r'$\omega_p t$', fontsize=14) 
    ax.set_ylabel(r'$k_p\sigma_x$', fontsize=14)
    if magn_check(np.sqrt(xsq)):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15) 
    
    saveas_eps_pdf(fig_sx, savepath, 'sigma_x_proj')                  
    plt.close(fig_sx)

    fig_sp = plt.figure()    
    plt.plot(t, np.sqrt(psq))
    ax = plt.gca()
    ax.set_xlabel(r'$\omega_p t$', fontsize=14) 
    ax.set_ylabel(r'$\sigma_{p_x}/m_e c$', fontsize=14)
    if magn_check(np.sqrt(psq)):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15)  
    
    saveas_eps_pdf(fig_sp, savepath, 'sigma_px_proj')                  
    plt.close(fig_sp)


    fig_xp = plt.figure()    
    plt.plot(t, xp)
    ax = plt.gca()
    ax.set_xlabel(r'$\omega_p t$', fontsize=14) 
    ax.set_ylabel(r'$ k_p \left \langle x\,p_x \right\rangle/m_e c$', fontsize=14)
    if magn_check(xp):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15)
    
    saveas_eps_pdf(fig_xp, savepath, 'xpx_proj')   
    plt.close(fig_xp)


    emittance_proj = np.sqrt(xsq*psq-np.power(xp,2))
    fig_e = plt.figure()    
    plt.plot(t, emittance_proj, label='projected emittance') 
    ax = plt.gca()
    ax.set_xlabel(r'$\omega_p t$', fontsize=14) 
    ax.set_ylabel(r'$k_p \epsilon_x$', fontsize=14)
    if magn_check(emittance_proj):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15)  
    
    saveas_eps_pdf(fig_e, savepath, 'emittance_proj')   
    plt.close(fig_e)

    fig_esl = plt.figure()    
    plt.plot(t, emittance_sliced, label='sliced emittance') 
    ax = plt.gca()
    ax.set_xlabel(r'$\omega_p t$', fontsize=14) 
    ax.set_ylabel(r'$k_p \epsilon_{x,\mathrm{sliced}}$', fontsize=14)
    if magn_check(emittance_sliced):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15)  
    
    saveas_eps_pdf(fig_esl, savepath, 'emittance_sliced')   
    plt.close(fig_esl)

    fig_e_slpr = plt.figure()    
    epp = plt.plot(t, emittance_proj, label='projected')
    esp = plt.plot(t, emittance_sliced, color = epp[0].get_color(), linestyle='--', label='sliced')    
    ax = plt.gca()
    ax.set_xlabel(r'$\omega_p t$', fontsize=14) 
    ax.set_ylabel(r'$k_p \epsilon_x$', fontsize=14)
    ax.legend()
    if magn_check(emittance_proj):
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15)  
    
    saveas_eps_pdf(fig_e_slpr, savepath, 'emittance_sl_proj')   
    plt.close(fig_e_slpr)


def plot_save_slice_rms_lines(slm, savepath, time = None, axdir=2, h5plot=True):

    if time == None:
        tidx = [0,-1];
    else:
        tidx = [(np.abs(slm.time_array - time)).argmin()]

    for i in tidx:

        if axdir == 2:
            sigma_xy = np.sqrt(slm.avgx2sq[i,:])
            sigma_xy_lab = r'$k_p \sigma_{x}$'
            sigma_xy_savename = 'sigma_x'
            sigma_pxy = np.sqrt(slm.avgp2sq[i,:])
            sigma_pxy_lab = r'$\sigma_{p_x}/mc$'
            sigma_pxy_savename = 'sigma_px'
            # Also define labels and savenames here!    
        elif axdir == 3:
            sigma_xy = np.sqrt(slm.avgx3sq[i,:])
            sigma_xy_lab = r'$k_p \sigma_{y}$'
            sigma_xy_savename = 'sigma_y'
            sigma_pxy = np.sqrt(slm.avgp3sq[i,:])
            sigma_pxy_lab = r'$\sigma_{p_y}/mc$'
            sigma_pxy_savename = 'sigma_py'

        fig_sigma_xy = plt.figure()
        plt.plot( slm.zeta_array,
                  sigma_xy)
        ax = plt.gca()
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
        ax.set_ylabel(sigma_xy_lab, fontsize=14)     
        saveas_eps_pdf(fig_sigma_xy, savepath, ('%s_time_%0.f' % (sigma_xy_savename, slm.time_array[i])))
        plt.close(fig_sigma_xy)


        fig_sigma_pxy = plt.figure()
        plt.plot( slm.zeta_array,
                  sigma_pxy)
        ax = plt.gca()
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
        ax.set_ylabel(sigma_xy_lab, fontsize=14)     
        saveas_eps_pdf(fig_sigma_pxy, savepath, ('%s_time_%0.f' % (sigma_pxy_savename, slm.time_array[i])))
        plt.close(fig_sigma_pxy)


def plot_save_slice_centroids(slm, savepath, h5plot=True):

    Xb0 = np.ones(slm.avgx2[0,:].shape)
    Yb0 = np.ones(slm.avgx3[0,:].shape)
    for i in range(0,len(slm.zeta_array)):
        if (slm.avgx2[0,i] != 0.0):
            Xb0[i] = slm.avgx2[0,i]
            Yb0[i] = slm.avgx3[0,i]


    Xb_norm = np.zeros( slm.avgx2.shape )
    Yb_norm = np.zeros( slm.avgx3.shape )
#    zeta_hseed = np.min(slm.zeta_array)
#    idx_hseed = (np.abs(slm.zeta_array-zeta_hseed)).argmin()

    for i in range(0,len(slm.zeta_array)):
#        if (slm.zeta_array[i] <= zeta_hseed):
        Xb_norm[:,i] = np.absolute( slm.avgx2[:,i]/Xb0[i] )
        Yb_norm[:,i] = np.absolute( slm.avgx3[:,i]/Yb0[i] )

    figXb = plt.figure()
    cax = plt.pcolormesh(   slm.zeta_array,
                            slm.time_array,
                            slm.avgx2,
                            cmap=cm.PuOr,
                            vmin=-np.amax(abs(slm.avgx2)), 
                            vmax=np.amax(abs(slm.avgx2)))
    cbar = figXb.colorbar(cax)   
    cbar.ax.set_ylabel('$k_p X_b$', fontsize=14)
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$\omega_p t$', fontsize=14) 
    saveas_png(figXb, savepath, 'Xb')
    plt.close(figXb)

    figYb = plt.figure()
    cax = plt.pcolormesh(   slm.zeta_array,
                            slm.time_array,
                            slm.avgx3,
                            cmap=cm.PuOr,
                            vmin=-np.amax(abs(slm.avgx3)), 
                            vmax=np.amax(abs(slm.avgx3)))
    cbar = figYb.colorbar(cax)
    cbar.ax.set_ylabel('$k_p Y_b$', fontsize=14)
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$\omega_p t$', fontsize=14)     
    saveas_png(figYb, savepath, 'Yb')
    plt.close(figYb)

    figXbnorm = plt.figure()
    cax = plt.pcolormesh(   slm.zeta_array,
                            slm.time_array,
                            Xb_norm,
                            cmap=cm.Blues,
                            vmin=0, 
                            vmax=np.amax(Xb_norm))
    cbar = figXbnorm.colorbar(cax)
    cbar.ax.set_ylabel('$|X_b/X_{b,0}|$', fontsize=14)
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$\omega_p t$', fontsize=14)  
    saveas_png(figXbnorm, savepath, 'Xb_rel')
    plt.close(figXbnorm)

    figYbnorm = plt.figure()
    cax = plt.pcolormesh(   slm.zeta_array,
                            slm.time_array,
                            Yb_norm,
                            cmap=cm.Blues,
                            vmin=0, 
                            vmax=np.amax(Yb_norm))
    cbar = figYbnorm.colorbar(cax)
    cbar.ax.set_ylabel('$|Y_b/Y_{b,0}|$', fontsize=14)
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$\omega_p t$', fontsize=14)
    saveas_png(figYbnorm, savepath, 'Yb_rel')
    plt.close(figYbnorm)

    figXb0 = plt.figure()
    plt.plot(slm.zeta_array, Xb0)
    ax = plt.gca()
    ymin, ymax = ax.get_ylim()
    if ymin > 0 and ymax > 0:
        plt.ylim(0, ymax*1.2)
    elif ymin < 0 and ymax < 0:
        plt.ylim(ymin*1.2, 0)        
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$k_p X_{b,0}$', fontsize=14)
    saveas_eps_pdf(figXb0, savepath, 'Xb0')
    plt.close(figXb0)

    figYb0 = plt.figure()
    plt.plot(slm.zeta_array, Yb0)
    ax = plt.gca()
    ymin, ymax = ax.get_ylim()
    if ymin > 0 and ymax > 0:
        plt.ylim(0, ymax*1.2)
    elif ymin < 0 and ymax < 0:
        plt.ylim(ymin*1.2, 0)    
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$k_p Y_{b,0}$', fontsize=14)     
    saveas_eps_pdf(figYb0, savepath, 'Yb0')
    plt.close(figYb0)

    figXbtail = plt.figure()
    plt.plot(slm.time_array, slm.avgx2[:,0])
    ax = plt.gca()
    ax.set_xlabel(r'$\omega_p t$', fontsize=14)
    ax.set_ylabel(r'$X_{b,\mathrm{tail}}$', fontsize=14)
    if magn_check(slm.avgx2[:,0]):    
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.gcf().subplots_adjust(left=0.18)
    else:
        plt.gcf().subplots_adjust(left=0.15)
    saveas_eps_pdf(figXbtail, savepath, 'Xb_tail')    
    plt.close(figXbtail)

    figYbtail = plt.figure()
    plt.plot(slm.time_array, slm.avgx3[:,0])
    ax = plt.gca()
    ax.set_xlabel(r'$\omega_p t$', fontsize=14)
    ax.set_ylabel(r'$Y_{b,\mathrm{tail}}$', fontsize=14) 
    saveas_eps_pdf(figYbtail, savepath, 'Yb_tail')   
    plt.close(figYbtail)  


def plot_save_slice_ene(slm, savepath):
    gamma = np.sqrt( 1 + np.power(slm.avgp1,2) 
                       + np.power(slm.avgp2,2)
                       + np.power(slm.avgp3,2) )

    fig = plt.figure()
    cax = plt.pcolormesh( slm.zeta_array,
                          slm.time_array,
                          gamma,
                          cmap=cm.GnBu,
                          vmin=0, vmax=np.amax(gamma))
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$\omega_p t$', fontsize=14)    
    cbar = fig.colorbar( cax )
    cbar.ax.set_ylabel(r'$\gamma$', fontsize=14)
    saveas_png(fig, savepath, 'gamma')
    plt.close(fig)    

def plot_curr_profile(slm, savepath):
    fig = plt.figure()
    dzeta = abs(slm.zeta_array[1] - slm.zeta_array[0]);
    curr = slm.charge[0,:] / dzeta

    # print('Q = %0.3e' % (np.sum(curr) * dzeta ))

    # I_A = 4 * pi * epsilon_0 * m * c^3 / e
    # [curr] = epsilon_0 * m * c^3 / e
    # curr * 4 * pi = I_b/I_A
    plt.plot(slm.zeta_array, curr/(4 * math.pi))
    ax = plt.gca()
    ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
    ax.set_ylabel(r'$I_{b,0}/I_A$', fontsize=14)
    saveas_eps_pdf(fig, savepath, 'Ib0')
    plt.close(fig)    

def main():

    parser = ps_parseopts()

    args = parser.parse_args()

    file = args.path

    slm = SliceMoms()
    slm.read(file)

    if args.zeta_range != None:
        slm.truncate_zeta_region(args.zeta_range[0], args.zeta_range[1])

    mkdirs_if_nexist(args.savepath)

    plot_curr_profile(slm, args.savepath)

    plot_save_slice_rms(slm, args.savepath)

    plot_save_slice_ene(slm, args.savepath)

    plot_save_proj_rms(slm, args.savepath, h5plot=args.h5plot)

    plot_save_slice_centroids(slm, args.savepath, h5plot=args.h5plot)

    plot_save_slice_rms_lines(slm, args.savepath, time = args.time, axdir=2, h5plot=args.h5plot)

if __name__ == "__main__":
    main()
