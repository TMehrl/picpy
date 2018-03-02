#!/usr/bin/env python3
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def main():
    ## QUANTITY_Slab_postrocessing    
    path='/Users/timonmehrling/ownCloud/PIC-codes/HiPACE/hipace/test/testSubGrid/DATA/'

    NHC=2
    if_incl_halo=1
    xoffset = 0.5

    file='mpsi_slab_0_nt_000000_proc_0'

    #cblims = [-1e-1 1e-1];
    
    summaryfile='output_summary.txt'

    grid_info = np.zeros((5, 3))
    
    with open(path + summaryfile, 'r') as csvfile:
        csvread = csv.reader(csvfile, delimiter='\t')
        i=0
        for line in csvread:
            for j in range(len(line)):
                grid_info[i,j] = float(line[j])
            i += 1        

    print(grid_info)

    Nx_wo_halo = int(grid_info[1,2])
    Ny_wo_halo = int(grid_info[2,2])
    
    Nx_incl_halo = Nx_wo_halo + 2*NHC;
    Ny_incl_halo = Ny_wo_halo + 2*NHC;

    dx=(grid_info[1,1] - grid_info[1,0]) / (Nx_wo_halo - 1)
    dy=(grid_info[2,1] - grid_info[2,0]) / (Ny_wo_halo - 1)
    
    x_array_wo_halo=np.linspace(grid_info[1,0],grid_info[1,1],Nx_wo_halo)
    y_array_wo_halo=np.linspace(grid_info[2,0],grid_info[2,1],Ny_wo_halo)
    
    x_array_incl_halo = np.linspace(grid_info[1,0]-2*dx,grid_info[1,1]+2*dx,Nx_incl_halo)
    y_array_incl_halo = np.linspace(grid_info[2,0]-2*dy,grid_info[2,1]+2*dy,Ny_incl_halo)
    
    filepath = path + file

    if if_incl_halo == 1:
        Nx=Nx_incl_halo
        Ny=Ny_incl_halo

        x_array=x_array_incl_halo
        y_array=y_array_incl_halo

    else:
        Nx=Nx_wo_halo
        Ny=Ny_wo_halo

        x_array=x_array_wo_halo
        y_array=y_array_wo_halo
    

    M_1D = np.fromfile(filepath,dtype=np.float32)
    M = np.transpose(M_1D.reshape((Nx, Ny)))

    fig = plt.figure()
    cax = plt.pcolormesh( x_array,
                          y_array,
                          M,
                          cmap=cm.PuBu)

    ax = plt.gca()
    ax.set_ylabel(r'$k_p y$', fontsize=14)
    ax.set_xlabel(r'$k_p x$', fontsize=14)
    cbar = fig.colorbar(cax)
    #cbar.ax.set_ylabel( gen_pretty_grid_name( self.g3d.name ), fontsize=14 )
    plt.show()    

    X,Y=np.meshgrid(x_array-xoffset,y_array)
    diff=M + (X ** 2 + Y ** 2) / 4
    fig = plt.figure()
    cax = plt.pcolormesh( x_array,
                          y_array,
                          diff,
                          cmap=cm.PuBu)

    ax = plt.gca()
    ax.set_ylabel(r'$k_p y$', fontsize=14)
    ax.set_xlabel(r'$k_p x$', fontsize=14)
    cbar = fig.colorbar(cax)
    plt.show()
    

if __name__ == "__main__":
    main()  