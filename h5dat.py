#!/usr/bin/env python3
# h5dat.py

import os
import sys
import numpy as np
import h5py
import picdefs


# Dictionary for PIC codes
piccodes = { 'hipace':'hipace',
              'osiris':'osiris'
            }


# General HDF5 file class with routines
# to check whether file is an HDF5 file
# and to print keys of attributes and datasets
class H5File:
    def __init__(self, file, h5ftype=None):

        # Allowed hdf5 extensions:
        self.__h5exts = ['.h5','.hdf5']
        # Grid 3D types in filenames:
        self.__g3dtypes = ['density', 'field', 'current']
        # RAW types in filenames:
        self.__rawtypes = ['raw']

        self.file = file
        self.h5ftype = h5ftype

    # Returning boolean: if file extension is hdf5 extension
    def is_h5_file(self, fext=None):
        if (fext!=None):
            return any(fext == h5ext for h5ext in self.__h5exts)
        else:
            fname, fext = os.path.splitext(self.file)
            return any(fext == h5ext for h5ext in self.__h5exts)          


    # Returning boolean: if file name contains 'raw'
    def is_g3d_file(self, fname):
        return any((mq in fname) for mq in self.__g3dtypes)

    # Returning boolean: if file name contains 'raw'
    def is_raw_file(self, fname):
        return any((mq in fname) for mq in self.__rawtypes)

    # Returning boolean:  if file extension is hdf5 extension and
    #                     if file name contains name of grid quantity
    def is_h5g3d_file(self):
        fname, fext = os.path.splitext(self.file)
        return self.is_h5_file(fext=fext) and self.is_g3d_file(fname=fname)

    # Returning boolean:  if file extension is hdf5 extension and
    #                     if file name contains 'raw'
    def is_h5raw_file(self):
        fname, fext = os.path.splitext(self.file)
        return self.is_h5_file(fext=fext) and self.is_raw_file(fname=fname)

    def fcheck(self):
        if self.h5ftype == 'g3d':
            return self.is_h5g3d_file()
        elif self.h5ftype == 'raw':       
            return self.is_h5raw_file()
        else:
            print('Error: No file type specified ["g3d", "raw"]!')
            sys.exit(1)

    def print_datasets(self):
        with h5py.File(self.file,'r') as hf:
            # Printing attributes
            print('HDF5 file datasets:')
            for item in hf.keys():
                print('\t' + item + ":", hf[item])

    def print_attributes(self):
        with h5py.File(self.file,'r') as hf:
            # Printing attributes
            print('HDF5 file attributes:')
            for item in hf.attrs.keys():
                print('\t' + item + ":", hf.attrs[item])

    def get_allowed_h5exts(self):
        return self.__h5exts

# Keys for HDF5 files
# Keys for PIC HDF5 files
class H5Keys:
    def __init__(self, piccode):

        # OSIRIS
        if piccode == piccodes['osiris']:

            # HDF5 GRID dataset keys
            self.__g3dkeys =  { 'density' : 'density'
                              }

            # HDF5 RAW dataset keys
            self.__rawkeys =  { 'x1':'x1',
                                'x2':'x2',
                                'x3':'x3',
                                'q':'q',
                                'p1':'p1',
                                'p2':'p2',
                                'p3':'p3'
                              }

            # HDF5 Attribute Keys
            self.__attrkeys = { 'nx':'NX',
                                'xmin':'XMIN',
                                'xmax':'XMAX',
                                'time':'TIME',
                                'dt':'DT',
                                'type':'TYPE',
                                'name':'NAME',
                              }
        # HiPACE
        elif piccode == piccodes['hipace']:

            # HDF5 GRID dataset keys
            self.__g3dkeys =  { 'beam_charge' : 'beam_charge',
                                'plasma_charge' : 'plasma_charge'
                              }
            # HDF5 GRID dataset types
            # Move these somewhere else...
            self.__g3dtypes = { 'density' : 'density',
                                'field' : 'field',
                                'current' : 'current'
                              }

            # HDF5 RAW dataset keys
            self.__rawkeys =  { 'x1':'x1',
                                'x2':'x2',
                                'x3':'x3',
                                'q':'q',
                                'p1':'p1',
                                'p2':'p2',
                                'p3':'p3'
                              }

            # HDF5 Attribute Keys
            self.__attrkeys = { 'nx':'NX',
                                'xmin':'XMIN',
                                'xmax':'XMAX',
                                'time':'TIME',
                                'dt':'DT',
                                'type':'TYPE',
                                'name':'NAME',
                              }


    def get_g3dkey(self,key):
        return self.__g3dkeys[key]
    def get_g3dkeys(self):
        return self.__g3dkeys
    def print_g3dkeys(self):
        for key in self.__g3dkeys: print(key)

    def get_rawkey(self,key):
        return self.__rawkeys[key]
    def get_rawkeys(self):
        return self.__rawkeys
    def print_rawkeys(self):
        for key in self.__rawkeys: print(key)

    def get_attrkey(self, key):
        return self.__attrkeys[key]
    def get_attrkeys(self):
        return self.__attrkeys
    def print_attrkeys(self):
        for key in self.__attrkeys: print(key)



class HiFile(H5Keys, H5File):
    def __init__(self, file):
        H5Keys.__init__(self, 'hipace')
        H5File.__init__(self, file)
        self.file = file

    def read_attrs(self):
        # Reading attributes
        with h5py.File(self.file,'r') as hf:
            self.nx = hf.attrs[   self.get_attrkey('nx') ]
            self.xmin = hf.attrs[ self.get_attrkey('xmin') ]
            self.xmax = hf.attrs[ self.get_attrkey('xmax') ]
            self.time = hf.attrs[ self.get_attrkey('time') ]
            self.dt = hf.attrs[   self.get_attrkey('dt') ]
            type_bytes = hf.attrs[ self.get_attrkey('type') ]
            name_bytes = hf.attrs[ self.get_attrkey('name') ]
            self.type = type_bytes[0].decode('UTF-8')
            self.name = name_bytes[0].decode('UTF-8')
            # Make these private and write getter functions!

    def get_x_arr(self,dim):
        return np.linspace(self.xmin[dim],self.xmax[dim],self.nx[dim])

    def get_z_arr(self):
        return np.linspace(self.time+self.xmin[0],self.time+self.xmax[0],self.nx[0])

    def get_zeta_arr(self):
        return np.linspace(self.xmin[0],self.xmax[0],self.nx[0])

    def get_xi_arr(self):
        return np.linspace(-self.xmin[0],-self.xmax[0],self.nx[0])

    def get_nt(self):
        return round(self.time/self.dt)

    def get_nx(self,dim):
        return self.nx[dim]        


# class H5PIC(H5File):
#   def __init__(self, file, piccode):
#
#     H5File.__init__(self, file)
#
#     self.piccode = piccode
#
#     if self.piccode == piccodes['hipace']:
#       self.file = HiH5File()
#     #elif self.piccode == piccodes['osiris']:
#       #self.file = OsH5File()
#     else:
#       print('Error:\tPIC code "%s" is not supported!' % piccode)
#       print('\tAllowed PIC codes: ',list(piccodes.values()))
#       sys.exit()


class HiRAW(HiFile):
    def __init__(self, file):
        HiFile.__init__(self, file)

    def read_data(self):

        with h5py.File(self.file,'r') as hf:

            # Reading datasets
            self.x1 = np.array(hf.get(  self.get_rawkey('x1') ))
            self.x2 = np.array(hf.get(  self.get_rawkey('x2') ))
            self.x3 = np.array(hf.get(  self.get_rawkey('x3') ))
            self.q =  np.array(hf.get(  self.get_rawkey('q')  ))
            self.p1 = np.array(hf.get(  self.get_rawkey('p1') ))
            self.p2 = np.array(hf.get(  self.get_rawkey('p2') ))
            self.p3 = np.array(hf.get(  self.get_rawkey('p3') ))

            self.npart = len(self.q)

#### 3D-grid
class Grid3d(HiFile):
    def __init__(self, file):
        HiFile.__init__(self, file)
        self.read_attrs()
        with h5py.File(self.file,'r') as hf:
            self.h5keys = list(hf.keys())
            if len(self.h5keys) == 1:
                self.dsetkey = self.h5keys[0]
            elif len(self.h5keys) == 0:
                print('Error:\tHDF5 file "%s" does not '
                  'contain any dataset!' %(self.file) )
                sys.exit()
            elif len(self.h5keys) > 1:
                print('Error:\tHDF5 file "%s" contains more '
                  'than one dataset!' %(self.file) )
                sys.exit()

        # self.fid = h5py.h5f.open(self.file.encode())
        # self.dset = h5py.h5d.open(self.fid, self.dsetkey.encode())

    def read_3D(self):
        with h5py.File(self.file,'r') as hf:
            # Reading dataset (here not caring how dataset is called)
            self.data = hf[self.dsetkey][()]


    def read_2D(self, i0=None, i1=None, i2=None):
        with h5py.File(self.file,'r') as hf:
        # Reading dataset (here not caring how dataset is called)
            if i0!=None and i1==None and i2==None:
                slice = hf[self.dsetkey][i0,:,:]
            elif i0==None and i1!=None and i2==None:
                slice = hf[self.dsetkey][:,i1,:]
            elif i0==None and i1==None and i2!=None:
                slice = hf[self.dsetkey][:,:,i2]
            else:
                print('Error:\tExactly one index must '
                  'be provided for HDF slice read in!')
                sys.exit()
        return(slice)

    def read_1D(self, i0=None, i1=None, i2=None):
        with h5py.File(self.file,'r') as hf:
            # Reading dataset (here not caring how dataset is called)
            if i0!=None and i1!=None and i2==None:
                line = hf[self.dsetkey][i0,i1,:]
            elif i0!=None and i1==None and i2!=None:
                line = hf[self.dsetkey][i0,:,i2]
            elif i0==None and i1!=None and i2!=None:
                line = hf[self.dsetkey][:,i1,i2]
            else:
                print('Error:\tExactly two indices must '
                  'be provided for HDF line read in!')
                sys.exit()
        return(line)

    def read_0D(self, i0=None, i1=None, i2=None):
        with h5py.File(self.file,'r') as hf:
            # Reading dataset (here not caring how dataset is called)
            if i0!=None and i1!=None and i2!=None:
                point = hf[self.dsetkey][i0,i1,i2]
            else:
                print('Error:\tExactly three indices must '
                  'be provided for HDF point read in!')
                sys.exit(1)
        return(point)

    def __n_none(self, arg0, arg1, arg2):
        return sum([arg0 == None, arg1 == None, arg2 == None])

    def read(self, 
             i0=None, 
             i1=None, 
             i2=None,
             x0=None,
             x1=None,
             x2=None):
        if self.__n_none(x0, x1, x2) == 3:
            if self.__n_none(i0, i1, i2) == 3:
                return self.read_3D()
            elif self.__n_none(i0, i1, i2) == 2:
                return self.read_2D(i0=i0, i1=i1, i2=i2) 
            elif self.__n_none(i0, i1, i2) == 1:
                return self.read_1D(i0=i0, i1=i1, i2=i2)             
            elif self.__n_none(i0, i1, i2) == 0:
                return self.read_0D(i0=i0, i1=i1, i2=i2)
        elif self.__n_none(i0, i1, i2) == 3:
            if x0 != None:
                i0 = np.argmin(np.abs(self.get_x_arr(0)-x0))
            if x1 != None:
                i1 = np.argmin(np.abs(self.get_x_arr(1)-x1))
            if x2 != None:
                i2 = np.argmin(np.abs(self.get_x_arr(2)-x2))
            return self.read(i0=i0, i1=i1, i2=i2)                    
        else:
            print('Error:\tMixed indices and positions '
              'not allowed for grid 3d read in!')
            sys.exit(1)                                  




class SliceMoms(H5File):
    def __init__(self, file):

        H5File.__init__(self, file)

        self.read()

    def read(self):
        with h5py.File(self.file,'r') as hf:
            # read order of moments
            order = 2
            # Reading datasets
            self.time_array = np.array(hf.get( 'time_array' ))
            self.zeta_array = np.array(hf.get( 'zeta_array' ))
            self.avgx1 = np.array(hf.get( 'avgx1' ))
            self.avgx2 = np.array(hf.get( 'avgx2' ))
            self.avgx3 = np.array(hf.get( 'avgx3' ))
            self.avgp1 = np.array(hf.get( 'avgp1' ))
            self.avgp2 = np.array(hf.get( 'avgp2' ))
            self.avgp3 = np.array(hf.get( 'avgp3' ))

            if order > 1:
                self.avgx1sq = np.array(hf.get( 'avgx1sq' ))
                self.avgx2sq = np.array(hf.get( 'avgx2sq' ))
                self.avgx3sq = np.array(hf.get( 'avgx3sq' ))
                self.avgp1sq = np.array(hf.get( 'avgp1sq' ))
                self.avgp2sq = np.array(hf.get( 'avgp2sq' ))
                self.avgp3sq = np.array(hf.get( 'avgp3sq' ))
                self.avgx1p1 = np.array(hf.get( 'avgx1p1' ))
                self.avgx2p2 = np.array(hf.get( 'avgx2p2' ))
                self.avgx3p3 = np.array(hf.get( 'avgx3p3' ))



class H5FList():
    def __init__(self, paths, h5ftype=None):

        self.paths = paths
        self.h5ftype = h5ftype

    def get(self, verbose=True):
        if not self.paths:
            print('Error: No file provided!')
            sys.exit(1)

        flist = []

        for path in self.paths:
            if os.path.isfile(path):
                file = path
                h5f = H5File(file, self.h5ftype)
                if h5f.fcheck():
                    flist.append(file)
                else:
                    if verbose: print('Skipping: ' + file)
            elif os.path.isdir(path):
                if verbose: print('"' + path + '"' + ' is a directory.')
                if verbose: 
                    print('Processing all ' + self.h5ftype + 
                          ' files in the provided directory.')
                for root, dirs, files in os.walk(path):
                    for filename in files:
                        file = root + '/' + filename
                        h5f = H5File(file, self.h5ftype)
                        if h5f.fcheck():
                            flist.append(file)
                        else:
                           if verbose: print('Skipping: ' + file)
            elif not os.path.exists(path):
                print('Error: Provided path does not exist!')
                sys.exit()
            else:
                print('Error: Provided path is neither a file nor a directory!')
                sys.exit()
        # Alphabetically sorting list
        flist = sorted(flist)
        return flist


def mkdirs_if_nexist( path ):
    folders = []

    while 1:
        path, folder = os.path.split(path)
        if folder != "":
            folders.append(folder)
        else:
            if path != "":
                folders.append(path)
            break
    folders.reverse()

    path = ""
    for folder in folders:
        path = path + folder + "/"
        if not os.path.isdir(path):
            print("Creating folder: " + path)
            os.mkdir(path)
