#!/usr/bin/env python3
# pp_h5dat.py

import os
import sys
import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt
from pp_h5dat import H5Plot
from pp_h5dat import mkdirs_if_nexist

def saveas_eps_pdf(fig, savepath, savename, h5plot=True, verbose=True, fformat='pdf'):

    fformat_list = ['eps','pdf']

    if fformat not in fformat_list:
        print("Error: fformat must be one of: " + fformat_list)

    mkdirs_if_nexist(savepath)
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


def saveas_png(fig, savepath, savename, verbose=True, dpi=600):

    fformat = 'png'
    mkdirs_if_nexist(savepath)
    fig.savefig( savepath + '/' + savename + '.' + fformat,
                  format=fformat,
                  dpi=dpi)    
      
    if verbose: 
        sys.stdout.write('Saved "%s.%s" at: %s/\n' % (savename, fformat, savepath))
        sys.stdout.flush() 

