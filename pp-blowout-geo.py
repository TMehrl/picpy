#!/usr/bin/env python3

import sys
import math
import argparse
import numpy as np
import matplotlib
import itertools
# Force matplotlib to not use any Xwindows backend.
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter
import scipy.constants as constants
import mpmath
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


class Blowout():
    def __init__(self, g3d, args):
      self.g3d = g3d
      self.args = args
      self.get_set_r_z_plane()


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
                             self.r_z_slice)

        ax = plt.gca()
        ax.set_ylabel(r'$k_p r$', fontsize=14)
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
        cbar = fig.colorbar(cax)

        filesuffix = '_t_%06.f' % (np.floor(self.g3d.time))
        fileprefix = 'plasma_charge'
        saveformat = 'png'
        savepath = self.args.savepath

        mkdirs_if_nexist(savepath)

        savename = fileprefix + filesuffix + '.' + saveformat

        fig.savefig( savepath + '/' + savename,
          format=saveformat,
          dpi=600)

        if self.args.verbose: print('Saved "' + savename + '" at: ' + self.args.savepath)


    def plot_r_lines(self):

        if not 'self.r_z_slice_model' in locals():
          self.gen_model_sheath()
            
        z_inds = np.int16(np.linspace(0,self.Nzeta-1,10))

        fig = plt.figure()
        for z_ind in z_inds:
            p = plt.plot(self.r_array, self.r_z_slice[:,z_ind])
            p[0].get_color()
            plt.plot(self.r_array, self.r_z_slice_model[:,z_ind],color = p[0].get_color(),linestyle='dashed')
            

        ax = plt.gca()
        ax.set_ylabel(r'$\rho$', fontsize=14)
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)

        filesuffix = '_t_%06.f' % (np.floor(self.g3d.time))
        fileprefix = 'plasma_charge_r'
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



    def calc_rb(self):

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
            
        slope = moving_average(np.asarray(slope), 50)

        fig = plt.figure()
        plt.plot(self.zeta_array,slope)
        ax = plt.gca()
        ax.set_ylabel('slope', fontsize=14)
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)           
        plt.show()

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


        filesuffix = '_t_%06.f' % (np.floor(self.g3d.time))
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

    def calc_deltarho_rhomax_sheathcharge(self):

        def f_delta_rho(Q,rho_peak,rb):
            if Q > 0 and rho_peak > 0:
                delta_rho = 0.5*(np.sqrt(4.0*Q/rho_peak + rb*rb) - rb)
            else:
                delta_rho = 0.01;
            return delta_rho

        idx_max = np.argmax(self.r_z_slice, axis=0)

        charge = np.zeros(self.Nzeta,dtype=np.float32)
        rho_max = np.zeros(self.Nzeta,dtype=np.float32)
        delta_rho = np.zeros(self.Nzeta,dtype=np.float32)

        for i in range(0,self.Nzeta):
            temp_array = self.r_array - self.rb[i]
            temp_array[temp_array<0] = np.amax(self.r_array)
            idx = (temp_array).argmin()
            charge[i] = np.sum(self.r_z_slice[idx:-1,i],axis=0)
            rho_max[i] = self.r_z_slice[idx_max[i],i]
            delta_rho[i] = f_delta_rho(Q=charge[i],rho_peak=rho_max[i],rb=self.rb[i])

        self.rho_max = moving_average(rho_max,4)
        self.charge = charge
        ########################## JUST TESTING:
        self.delta_rho = moving_average(delta_rho, 4) * 0.1
        ########################## JUST TESTING, need to improve estimate for delta_rho!!

    def plot_deltarho_rhomax_sheathcharge(self):

        if not 'self.rb' in locals():
          self.calc_rb()

        if not 'self.rho_max' in locals():
          self.calc_deltarho_rhomax_sheathcharge()

        fig_charge = plt.figure()
        plt.plot(self.zeta_array, self.charge)
        ax_charge = plt.gca()
        ax_charge.set_ylabel(r'$Q$', fontsize=14)
        ax_charge.set_xlabel(r'$k_p \zeta$', fontsize=14)    
        charge_fileprefix = 'charge'

        fig_rhomax= plt.figure()
        plt.plot(self.zeta_array, self.rho_max)
        ax_rhomax = plt.gca()
        ax_rhomax.set_ylabel(r'$\rho_\mathrm{max}$', fontsize=14)
        ax_rhomax.set_xlabel(r'$k_p \zeta$', fontsize=14) 
        rhomax_fileprefix = 'rho_max'

        fig_deltarho = plt.figure()
        plt.plot(self.zeta_array, self.delta_rho)
        ax_deltarho = plt.gca()
        ax_deltarho.set_ylabel(r'$\Delta_\rho$', fontsize=14)
        ax_deltarho.set_xlabel(r'$k_p \zeta$', fontsize=14) 
        deltarho_fileprefix = 'delta_rho'

        filesuffix = '_t_%06.f' % (np.floor(self.g3d.time))
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

        if not 'self.rb' in locals():
          self.calc_rb()

        if not 'self.rho_max' in locals():
          self.calc_deltarho_rhomax_sheathcharge()

        self.r_z_slice_model = -1*np.ones(self.r_z_slice.shape,dtype=np.float32)
        for i in range(0,self.Nzeta):
            ridx = np.where(self.r_array > self.rb[i])
            self.r_z_slice_model[ridx,i] = f_rho(self.r_array[ridx],self.rb[i],self.rho_max[i],self.delta_rho[i])

    def plot_model_sheath(self):

        if not 'self.r_z_slice_model' in locals():
          self.gen_model_sheath()
            
        fig = plt.figure() 
        cax = plt.pcolormesh(self.zeta_array,
                             self.r_array,
                             self.r_z_slice_model)

        ax = plt.gca()
        ax.set_ylabel(r'$k_p r$', fontsize=14)
        ax.set_xlabel(r'$k_p \zeta$', fontsize=14)
        cbar = fig.colorbar(cax)

        filesuffix = '_t_%06.f' % (np.floor(self.g3d.time))
        fileprefix = 'plasma_charge_model'
        saveformat = 'png'
        savepath = self.args.savepath

        mkdirs_if_nexist(savepath)

        savename = fileprefix + filesuffix + '.' + saveformat

        fig.savefig( savepath + '/' + savename,
          format=saveformat,
          dpi=600)

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
