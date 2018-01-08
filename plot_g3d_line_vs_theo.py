#!/usr/bin/env python3

import math
import numpy as np
import matplotlib.pyplot as plt
import plot_g3d

def sim_theo_cmp_plot( g3d_p ):
    y_theo = np.cos(g3d_p.x_array)

    fig = plt.figure()
    cax = plt.plot( g3d_p.x_array, y_theo)
    cax = plt.plot( g3d_p.x_array, g3d_p.line)
    plt.show()

def main():
    parser = plot_g3d.parser( ptype='line' )
    args = parser.parse_args()
    
    flist = plot_g3d.gen_filelist( args )

    for file in flist:
        g3d_p = plot_g3d.G3d_plot_line(file, args)
        sim_theo_cmp_plot(g3d_p)

if __name__ == "__main__":
    main()
