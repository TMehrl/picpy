#!/usr/bin/env python3

import csv
import os
import glob
import sys
import numpy as np
import argparse

import h5py

import matplotlib.pyplot as plt
from matplotlib import cm
from pp_h5dat import mkdirs_if_nexist

def binSlab_parser():
    
    desc = """This is the picpy postprocessing tool."""

    parser = argparse.ArgumentParser( description=desc)
    parser.add_argument(  'path',
                          metavar = 'PATH',
                          nargs = '*',
                          help = 'Path to binary files.')
    parser.add_argument(  "-o", "--output-path",
                          action="store",
                          dest="output_path",
                          metavar="PATH",
                          default=False,
                          help = """Path to which generated files will be saved.
                              (default: %(default)s)""")
    parser.add_argument(  "-k", "--keep",
                          action="store_true",
                          dest="keep",
                          default=False,
                          help = """Keep the old binary files 
                                instead of deleting them""")
    return parser
    
    
def main():
    
    parser = binSlab_parser() 
    args = parser.parse_args()
    


    
    for files in args.path:
        
        filename = ''
        slash = '/'
        if slash in files:
            dirname, filename = files.split('/')
        else:
            dirname = filepath
            
        array = np.empty(0)
        N = np.empty(0)
        box_dimensions = np.empty(0)
        data_type = 'none'
        output_name = ''
        name = ''
        if args.output_path:
            output_path = args.output_path
        else: 
            output_path = dirname

        timestamp = files.split("_")[-1]
        
        if 'field' in files:
            output_name += 'field'
            data_type = 'field'
            if 'Ez' in files:
                output_name += '_Ez'
                name = 'Ez'
                print('its the ez field!')
            elif 'ExmBy' in files:
                output_name += '_ExmBy'
                name = 'ExmBy'
                print('its the ExmBy field!')
                
            
        else:
            data_type = 'density'
            output_name += 'density_plasma_electrons'
            name = 'density_plasma_electrons'
            print('its a density!')
    
        array = np.fromfile(files, dtype=np.float64)
        Nz = np.int(array[0])
        Nr = np.int(array[1])
        # print('Nz: %i , Nr: %i' %(Nz, Nr))
        
        zmin, rmin, zmax, rmax = array[2:6]
        
        data = np.reshape(array[6:], (Nz, Nr))
        #data = np.transpose(data)
        # print(np.shape(data))
        if 'ExmBy' in files:
            data = np.append(-np.flip(data[:,1:], axis=1), data, axis=1)
        else:
            data = np.append(np.flip(data[:,1:], axis=1), data, axis=1)
        # print(np.shape(data))
        # control plot if needed
        # plt.pcolormesh(data, cmap=cm.plasma)
        # plt.show()
        print(output_path + output_name + '_' + timestamp + '.h5')
        mkdirs_if_nexist(output_path)
        f = h5py.File(output_path + '/' + output_name + '_' + timestamp + '.h5', 'w')
        dset = f.create_dataset( output_name + '_' + timestamp, shape=np.shape(data), data=data)
        f.attrs['DT'] = [np.float32(0.0)]
        #f.create_dataset('NAME',(32,), dtype="S10")
        f.attrs['NAME'] = name#repr(output_name) #.encode('UTF-8')
        #dt = h5py.special_dtype(vlen=str)     # PY3
        #dset = f.create_dataset('NAME', dtype=dt, data=output_name)
        f.attrs['NX'] = [np.int32(np.shape(data)[0]), np.int32(np.shape(data)[1]) ]
        f.attrs['TIME'] = [np.float32(timestamp)]
        f.attrs['TYPE'] = data_type #repr(data_type) #.encode('UTF-8')
        #dset = f.create_dataset('TYPE', dtype=dt, data=data_type)
        f.attrs['XMAX'] = [np.float32(zmax), np.float32(rmax)]
        f.attrs['XMIN'] = [np.float32(zmin), np.float32(-rmax)]
        print('Writing new hdf5 file...')
        f.flush()
        
        if not args.keep:
            print('Deleting old binary file...')
            #os.remove(files)
        print('Done!')
        
if __name__ == "__main__":
    main() 
    
    
    
    