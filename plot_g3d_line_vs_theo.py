#!/usr/bin/env python3

import math
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import plot_g3d

def sim_theo_cmp_plot( g3d_p ):
    y_theo = np.max(g3d_p.line) * np.sin(-1*g3d_p.x_array)

    fig = plt.figure()
    cax = plt.plot( g3d_p.x_array, g3d_p.line)
    cax = plt.plot( g3d_p.x_array, y_theo, '--')

    g3d_p.mkdirs_if_nexist()

    saveformat = g3d_p.args.file_format
    filesuffix = '_%06.f' % (np.floor(g3d_p.g3d.time))

    fileprefix = g3d_p.g3d.name

    savename = fileprefix + filesuffix + '.' + saveformat

    if saveformat=='png':
        fig.savefig(  g3d_p.args.savepath + '/' + savename,
                  format=saveformat,
                  dpi=600)
    else:
        fig.savefig(  g3d_p.args.savepath + '/' + savename,
                      format=saveformat)
    if g3d_p.args.verbose: print('Saved "' + savename + '" at: ' + g3d_p.args.savepath)

    if g3d_p.args.ifshow: plt.show()
    plt.close(fig)



def main():
    parser = plot_g3d.parser( ptype='line' )
    args = parser.parse_args()
    
    flist = plot_g3d.gen_filelist( args )

    for file in flist:
        g3d_p = plot_g3d.G3d_plot_line(file, args)
        sim_theo_cmp_plot(g3d_p)

if __name__ == "__main__":
    main()
