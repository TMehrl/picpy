#!/usr/bin/env python3
# pp_h5dat.py

import os
import sys
import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from pp_h5dat import H5Plot
from pp_h5dat import mkdirs_if_nexist

colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

# Sort colors by hue, saturation, value and name.
color_names_sorted = [name for hsv, name in 
                        sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name)
                                for name, color in colors.items())]


def show_color_palette():
    n = len(color_names_sorted)
    ncols = 4
    nrows = n // ncols + 1

    fig, ax = plt.subplots(figsize=(8, 5))

    # Get height and width
    X, Y = fig.get_dpi() * fig.get_size_inches()
    h = Y / (nrows + 1)
    w = X / ncols

    for i, name in enumerate(color_names_sorted):
        col = i % ncols
        row = i // ncols
        y = Y - (row * h) - h

        xi_line = w * (col + 0.05)
        xf_line = w * (col + 0.25)
        xi_text = w * (col + 0.3)

        ax.text(xi_text, y, name, fontsize=(h * 0.8),
                horizontalalignment='left',
                verticalalignment='center')

        ax.hlines(y + h * 0.1, xi_line, xf_line,
                  color=colors[name], linewidth=(h * 0.6))

    ax.set_xlim(0, X)
    ax.set_ylim(0, Y)
    ax.set_axis_off()

    fig.subplots_adjust(left=0, right=1,
                        top=1, bottom=0,
                        hspace=0, wspace=0)
    plt.show()


def saveas_eps_pdf(fig, savepath, savename, h5plot=True, verbose=True, fformat='pdf', transparent=True):

    fformat_list = ['eps','pdf']

    if fformat not in fformat_list:
        print("Error: fformat must be one of: " + fformat_list)

    mkdirs_if_nexist(savepath)
    fig.savefig( savepath + '/' + savename + '.' + fformat,
                  format=fformat, transparent=transparent)
      
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

