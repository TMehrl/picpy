#!/usr/bin/env python3

import sys
import math
import argparse
import numpy as np
import scipy.signal as ssignal
import matplotlib
import itertools
# Force matplotlib to not use any Xwindows backend.
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter
import scipy.constants as constants
import mpmath

mypath = os.path.dirname(os.path.abspath( __file__ ))
incpath = os.path.split(mypath)[0] + '/inc'
sys.path.append(incpath)
from pp_h5dat import Grid3d
from pp_h5dat import H5FList
from pp_h5dat import H5Plot
from pp_h5dat import mkdirs_if_nexist

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
        prefix = 'g3d_name'
        path = './plots'


def blowout_geo_parser():

    desc = """This is the picpy postprocessing tool."""
    # Line vs-theo plot arguments
    parser = argparse.ArgumentParser( description=desc)
    parser.add_argument(  'path',
                          metavar = 'PATH',
                          nargs = '*',
                          help = 'Path to grid file.')
    parser.add_argument(  "-v", "--verbose",
                          dest = "verbose",
                          action="store_true",
                          default=True,
                          help = "Print info (Default).")
    parser.add_argument(  "-q", "--quiet",
                          dest = "verbose",
                          action = "store_false",
                          help = "Don't print info.")    
    parser.add_argument(  "-s", "--save-path",
                          action="store",
                          dest="savepath",
                          metavar="PATH",
                          default=parsedefs.save.path + '/blowout-geo',
                          help = """Path to which generated files will be saved.
                              (Default: %(default)s)""")
    parser.add_argument(  "-f", "--format",
                          action='store',
                          dest="file_format",
                          metavar="FORMAT",
                          choices=[ 'png',
                                    'pdf',
                                    'eps',],
                          default='eps',
                          help= """Format of output file (Default: %(default)s).""")
    parser.add_argument(  "--h5",
                          dest = "h5plot",
                          action="store_true",
                          default=True,
                          help = "Save plot as hdf5 file (Default: %(default)s).")      
    return parser

def moving_average(data_set, periods=3):
    weights = np.ones(periods) / periods
    return np.convolve(data_set, weights, mode='same')

def f_rho(r_array,rb,rho_peak,delta_rho):
    return rho_peak * np.exp(-(r_array-rb)/delta_rho)


def floc_max(array, width):
    lm = ssignal.argrelmax(array, order=width)
    if lm[0].any():
        flm = lm[0][0]
    else:
        flm = 0
    return flm    


def fidx_loc_max(array, width):
    idx_lm = ssignal.find_peaks_cwt(array, np.arange(1,width))
    if idx_lm.size > 0:
        fidx_lm = idx_lm[0]
    else:
        fidx_lm = 0
    return fidx_lm  

class Blowout():
    def __init__(self, g3d, args):
      self.g3d = g3d
      self.args = args
      self.get_set_r_z_plane()
      self.mm_width = 10

    def get_set_r_z_plane(self):

        slice = -1*np.transpose(self.g3d.read(x2=0.0))
        x_array = self.g3d.get_x_arr(1)
        self.zeta_array = self.g3d.get_x_arr(0)
        self.Nzeta = self.g3d.get_nx(0)

        idx_gtr0 = x_array>=0.0

        self.r_array = x_array[idx_gtr0]
        self.Nr = len(idx_gtr0)
        self.r_z_slice = slice[idx_gtr0,:]
        self.dr = self.g3d.get_dx(1)

    def plot_r_z_plane(self):

        fig = plt.figure() 
        cax = plt.pcolormesh(self.zeta_array,
                             self.r_array,
                             self.r_z_slice,
                             cmap=plt.get_cmap('magma'))

        cbar = fig.colorbar(cax)
        cbar.ax.set_ylabel(r'$\rho$', fontsize=14 )

        ax = plt.gca()
        ax.set_ylabel(r'$k_p r$', fontsize=14)
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)

        filesuffix = '_t_%06.f' % (np.floor(self.g3d.get_time()))
        fileprefix = 'plasma_charge'
        saveformat = 'png'
        savepath = self.args.savepath

        mkdirs_if_nexist(savepath)

        savename = fileprefix + filesuffix + '.' + saveformat

        fig.savefig( savepath + '/' + savename,
          format=saveformat,
          dpi=600)
        plt.close(fig) 

        if self.args.verbose: print('Saved "' + savename + '" at: ' + self.args.savepath)


    def plot_r_lines(self):

        if not hasattr(self,'r_z_slice_model'):
          self.gen_model_sheath()
            
        z_inds = np.int16(np.linspace(0,self.Nzeta-1,10))

        fig = plt.figure()
        for z_ind in z_inds:
            p = plt.plot(self.r_array, self.r_z_slice[:,z_ind])
            plt.plot(self.r_array, self.r_z_slice_model[:,z_ind],color = p[0].get_color(),linestyle='--') 
            plt.plot(self.r_array, self.r_z_slice_model_flocmax[:,z_ind],color = p[0].get_color(),linestyle='-.')            

        ax = plt.gca()
        ax.set_ylabel(r'$\rho$', fontsize=14)
        ax.set_xlabel(r'$k_p r$', fontsize=14)

        filesuffix = '_t_%06.f' % (np.floor(self.g3d.get_time()))
        fileprefix = 'plasma_charge_r'
        saveformat = 'eps'
        savepath = self.args.savepath

        mkdirs_if_nexist(savepath)

        savename = fileprefix + filesuffix + '.' + saveformat

        fig.savefig(  savepath + '/' + savename,
                          format=saveformat)
        plt.close(fig) 

        if self.args.h5plot: 
            h5lp = H5Plot()
            h5lp.inherit_matplotlib_line_plots(ax)
            h5lp.write(savepath + '/' + fileprefix + filesuffix + '.h5')

        if self.args.verbose: print('Saved "' + savename + '" at: ' + self.args.savepath)



    def calc_rb(self):
        print('Calculating rb...')
        idx_max = np.argmax(self.r_z_slice, axis=0)

        rb = np.zeros(self.Nzeta,dtype=np.float32)
        slope = np.zeros(self.Nzeta,dtype=np.float32)

        for i in range(0,self.Nzeta):
            idx_gtr0_all = np.where(self.r_z_slice[:,i]>0)[0]
            if idx_gtr0_all.size != 0:
                if idx_gtr0_all[0] != 0:
                    idx_gtr0 = idx_gtr0_all[0]
                else:
                    idx_gtr0 = 1
            else:
                idx_gtr0 = 1

            if abs(self.r_z_slice[idx_gtr0,i]) > 0:
                ratio = abs(self.r_z_slice[idx_gtr0-1,i]) / abs(self.r_z_slice[idx_gtr0,i])
                slope[i] = (self.r_z_slice[idx_gtr0,i] - self.r_z_slice[idx_gtr0-1,i])/self.dr
            else:
                ratio = 0
                slope[i] = 0

            rb[i] = self.r_array[idx_gtr0-1] + self.dr * ratio /(1.0+ratio)
            
        slope = moving_average(np.asarray(slope), self.mm_width)

        fig = plt.figure()
        plt.plot(self.zeta_array,slope)
        ax = plt.gca()
        ax.set_ylabel('slope', fontsize=14)
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)           
        plt.show()
        plt.close(fig) 
        #self.rb = rb * (1.0 - np.exp(-0.5* np.power(slope/0.05,6))) 
        self.rb = rb 
        self.rb_maxrho = self.r_array[idx_max]       

    def plot_rb(self):

        fig = plt.figure()
        plt.plot(self.zeta_array,self.rb)
        plt.plot(self.zeta_array,self.rb_maxrho)

        ax = plt.gca()
        ax.set_ylabel(r'$k_p r_b$', fontsize=14)
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)    
        ax.legend([r'$r_b$ from zero crossing', r'$r_b$ from global maximum'])

        filesuffix = '_t_%06.f' % (np.floor(self.g3d.get_time()))
        fileprefix = 'plasma_charge_rb'
        saveformat = 'eps'
        savepath = self.args.savepath

        mkdirs_if_nexist(savepath)

        savename = fileprefix + filesuffix + '.' + saveformat

        fig.savefig(  savepath + '/' + savename,
                          format=saveformat)
        
        if self.args.h5plot: 
            h5lp = H5Plot()
            h5lp.inherit_matplotlib_line_plots(ax)
            h5lp.write(savepath + '/' + fileprefix + filesuffix + '.h5')

        if self.args.verbose: print('Saved "' + savename + '" at: ' + self.args.savepath)
        plt.close(fig) 

    def calc_deltarho_rhomax_sheathcharge(self,zeta_avg_region=[-3.0, 0.0]):

        def f_delta_rho(Q,rho_peak,rb):
            if Q > 0 and rho_peak > 0:
                delta_rho = 0.5*(np.sqrt(4.0*Q/rho_peak + rb*rb) - rb)
            else:
                delta_rho = 0.01;
            return delta_rho

        # Index of absolute maximum
        idx_max = np.argmax(self.r_z_slice, axis=0)

        charge = np.zeros(self.Nzeta,dtype=np.float32)
        rho_max = np.zeros(self.Nzeta,dtype=np.float32)
        rho_flocmax = np.zeros(self.Nzeta,dtype=np.float32)        
        delta_rho = np.zeros(self.Nzeta,dtype=np.float32)
        delta_rho_flocmax = np.zeros(self.Nzeta,dtype=np.float32)

        width = 5

        for i in range(0,self.Nzeta):
            temp_array = self.r_array - self.rb[i]
            temp_array[temp_array<0] = np.amax(self.r_array)
            idx = (temp_array).argmin()
            
            # Compute charge in sheath
            charge[i] = np.dot(self.r_z_slice[idx:-1,i], self.r_array[idx:-1]) *  self.dr
            rho_max[i] = self.r_z_slice[idx_max[i],i]
            delta_rho[i] = f_delta_rho(Q=charge[i],rho_peak=rho_max[i],rb=self.rb[i])
            
            # Indices of first local maxima
            fidx_lm = fidx_loc_max(self.r_z_slice[idx-width:-1,i], width)
            rho_flocmax[i] = self.r_z_slice[idx-width+fidx_lm,i]
            delta_rho_flocmax[i] = f_delta_rho(Q=charge[i],rho_peak=rho_flocmax[i],rb=self.rb[i])
            # rho_flocmax[i] = floc_max(self.r_z_slice[idx-width:-1,i], width)
            # delta_rho_flocmax[i] = f_delta_rho(Q=charge[i],rho_peak=rho_flocmax[i],rb=self.rb[i])

        self.rho_max = moving_average(rho_max, self.mm_width)
        self.rho_flocmax = moving_average(rho_flocmax, self.mm_width)
        self.charge = moving_average(charge, self.mm_width)
        self.delta_rho = moving_average(delta_rho, self.mm_width)
        self.delta_rho_flocmax = moving_average(delta_rho_flocmax, self.mm_width)

        if zeta_avg_region != [] and len(zeta_avg_region) == 2:
            region_idx = np.where((self.zeta_array >= zeta_avg_region[0]) & (self.zeta_array < zeta_avg_region[1]))
            avg_rho_max = np.mean(self.rho_max[region_idx])
            avg_rho_flocmax = np.mean(self.rho_flocmax[region_idx])
            avg_charge = np.mean(self.charge[region_idx])
            avg_delta_rho = np.mean(self.delta_rho[region_idx])
            avg_delta_rho_flocmax = np.mean(self.delta_rho_flocmax[region_idx])
            print('Zeta region: [%0.2f, %0.2f]' % (zeta_avg_region[0], zeta_avg_region[1]))
            print('Avg rho_max: %0.3f' % avg_rho_max) 
            print('Avg first local rho_max: %0.3f' % avg_rho_flocmax) 
            print('Avg charge: %0.3f' % avg_charge) 
            print('Avg Delta_rho: %0.3f' % avg_delta_rho)
            print('Avg Delta_rho using first local max: %0.3f' % avg_delta_rho_flocmax) 

    def plot_deltarho_rhomax_sheathcharge(self):

        if not hasattr(self,'rb'):
          self.calc_rb()

        if not hasattr(self,'rho_max'):            
          self.calc_deltarho_rhomax_sheathcharge()

        fig_charge = plt.figure()
        plt.plot(self.zeta_array, self.charge)
        ax_charge = plt.gca()
        ax_charge.set_ylabel(r'$\lambda=\int_{r_b}^\infty r \rho(r) dr$', fontsize=14)
        ax_charge.set_xlabel(r'$k_p \zeta$', fontsize=14)    
        charge_fileprefix = 'charge'

        fig_rhomax= plt.figure()
        plt.plot(self.zeta_array, self.rho_max)
        plt.plot(self.zeta_array, self.rho_flocmax)
        ax_rhomax = plt.gca()
        ax_rhomax.legend([r'using global maximum', r'using first local maximum'])
        ax_rhomax.set_ylabel(r'$\rho_\mathrm{max}$', fontsize=14)
        ax_rhomax.set_xlabel(r'$k_p \zeta$', fontsize=14) 
        rhomax_fileprefix = 'rho_max'

        fig_deltarho = plt.figure()
        plt.plot(self.zeta_array, self.delta_rho)
        plt.plot(self.zeta_array, self.delta_rho_flocmax)
        ax_deltarho = plt.gca()
        ax_deltarho.legend([r'using global maximum', r'using first local maximum'])
        ax_deltarho.set_ylabel(r'$\Delta_\rho$', fontsize=14)
        ax_deltarho.set_xlabel(r'$k_p \zeta$', fontsize=14) 
        deltarho_fileprefix = 'delta_rho'

        filesuffix = '_t_%06.f' % (np.floor(self.g3d.get_time()))
        saveformat = 'eps'
        savepath = self.args.savepath

        mkdirs_if_nexist(savepath)

        charge_savename = charge_fileprefix + filesuffix + '.' + saveformat
        rhomax_savename = rhomax_fileprefix + filesuffix + '.' + saveformat
        deltarho_savename = deltarho_fileprefix + filesuffix + '.' + saveformat

        fig_charge.savefig(savepath + '/' + charge_savename,
                           format=saveformat)
        if self.args.verbose: print('Saved "' + charge_savename + '" at: ' + self.args.savepath)   

        fig_rhomax.savefig(savepath + '/' + rhomax_savename,
                           format=saveformat)
        if self.args.verbose: print('Saved "' + rhomax_savename + '" at: ' + self.args.savepath)   

        fig_deltarho.savefig(savepath + '/' + deltarho_savename,
                           format=saveformat)
        if self.args.verbose: print('Saved "' + deltarho_savename + '" at: ' + self.args.savepath)   

        if self.args.h5plot: 
            h5lp = H5Plot()
            h5lp.inherit_matplotlib_line_plots(ax_charge)
            h5lp.write(savepath + '/' + charge_fileprefix + filesuffix + '.h5')

            h5lp = H5Plot()
            h5lp.inherit_matplotlib_line_plots(ax_rhomax)
            h5lp.write(savepath + '/' + rhomax_fileprefix + filesuffix + '.h5')

            h5lp = H5Plot()
            h5lp.inherit_matplotlib_line_plots(ax_deltarho)
            h5lp.write(savepath + '/' + deltarho_fileprefix + filesuffix + '.h5')

        plt.show()

    def gen_model_sheath(self):

        if not hasattr(self,'rb'):             
          self.calc_rb()

        if not hasattr(self,'rho_max'):
          self.calc_deltarho_rhomax_sheathcharge()

        self.r_z_slice_model = -1*np.ones(self.r_z_slice.shape,dtype=np.float32)
        self.r_z_slice_model_flocmax = -1*np.ones(self.r_z_slice.shape,dtype=np.float32)
        for i in range(0,self.Nzeta):
            ridx = np.where(self.r_array > self.rb[i])
            self.r_z_slice_model[ridx,i] = f_rho(self.r_array[ridx],self.rb[i],self.rho_max[i],self.delta_rho[i])
            self.r_z_slice_model_flocmax[ridx,i] = f_rho(self.r_array[ridx],self.rb[i],self.rho_flocmax[i],self.delta_rho_flocmax[i])

    def plot_model_sheath(self):

        if not hasattr(self,'r_z_slice_model'):          
          self.gen_model_sheath()
            
        fig = plt.figure() 
        cax = plt.pcolormesh(self.zeta_array,
                             self.r_array,
                             self.r_z_slice_model,
                             cmap=plt.get_cmap('magma'))

        ax = plt.gca()
        ax.set_ylabel(r'$k_p r$', fontsize=14)
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
        cbar = fig.colorbar(cax)
        cbar.ax.set_ylabel(r'$\rho$', fontsize=14 )

        filesuffix = '_t_%06.f' % (np.floor(self.g3d.get_time()))
        fileprefix = 'plasma_charge_model'
        saveformat = 'png'
        savepath = self.args.savepath

        mkdirs_if_nexist(savepath)

        savename = fileprefix + filesuffix + '.' + saveformat

        fig.savefig( savepath + '/' + savename,
          format=saveformat,
          dpi=600)
        plt.close(fig) 
        if self.args.verbose: print('Saved "' + savename + '" at: ' + self.args.savepath)


        fig = plt.figure() 
        cax = plt.pcolormesh(self.zeta_array,
                             self.r_array,
                             self.r_z_slice_model_flocmax,
                             cmap=plt.get_cmap('magma'))

        ax = plt.gca()
        ax.set_ylabel(r'$k_p r$', fontsize=14)
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
        cbar = fig.colorbar(cax)
        cbar.ax.set_ylabel(r'$\rho$', fontsize=14 )

        filesuffix = '_t_%06.f' % (np.floor(self.g3d.get_time()))
        fileprefix = 'plasma_charge_model_flocmax'
        saveformat = 'png'
        savepath = self.args.savepath

        savename = fileprefix + filesuffix + '.' + saveformat

        fig.savefig( savepath + '/' + savename,
          format=saveformat,
          dpi=600)
        plt.close(fig) 
        if self.args.verbose: print('Saved "' + savename + '" at: ' + self.args.savepath)

def main():

    parser = blowout_geo_parser()
    args = parser.parse_args()

    h5flist = H5FList(args.path, 'g3d')
    flist = h5flist.get()

    for file in flist:
        g3d = Grid3d(file)
        if (g3d.name == 'plasma_charge'):
            blowout = Blowout(g3d, args)
            blowout.plot_r_z_plane()
            blowout.plot_r_lines()
            blowout.plot_rb()
            blowout.calc_deltarho_rhomax_sheathcharge()
            blowout.plot_deltarho_rhomax_sheathcharge()
            blowout.plot_model_sheath()

        else:
            print('ERROR: HDF5 file must contain a plasma_charge dataset!')
            sys.exit(1)    


if __name__ == "__main__":
    main()
