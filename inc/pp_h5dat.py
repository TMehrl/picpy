#!/usr/bin/env python3
# pp_h5dat.py

import os
import sys
import math
import numpy as np
import re
import h5py


# Dictionary for PIC codes
piccodes = { 'hipace':'hipace',
              'osiris':'osiris'
            }


class File:
    def __init__(self, file, intent):    
        self.__file = file

        if 'r' in intent:
            if not self.is_file():
                print('Error:\tFile "%s" does not exist!' %(self.get_file()) )
                sys.exit(1)          

    def get_file(self):
        return self.__file

    def get_filename(self):
        return os.path.split(self.get_file())[1]  

    def get_path(self):
        return os.path.split(self.get_file())[0]

    def is_file(self):
        return os.path.isfile(self.get_file()) 



# General HDF5 file class with routines
# to check whether file is an HDF5 file
# and to print keys of attributes and datasets
class H5File(File):
    def __init__(self, file, intent):

        File.__init__(self, file, intent)

        # Allowed hdf5 extensions:
        self.__h5exts = ['.h5','.hdf5']

        if not self.is_h5_file():
            print('File "%s" is not an HDF5 file!' %(self.get_file()) )

    # Returning boolean: if file extension is hdf5 extension
    def is_h5_file(self, fext=None):
        if (fext!=None):
            return any(fext == h5ext for h5ext in self.__h5exts)
        else:
            fname, fext = os.path.splitext(self.get_file())
            return any(fext == h5ext for h5ext in self.__h5exts)          

    def print_datasets(self):
        with h5py.File(self.get_file(),'r') as hf:
            # Printing attributes
            print('HDF5 file datasets:')
            for item in hf.keys():
                print('\t' + item + ":", hf[item])

    def print_attributes(self):
        with h5py.File(self.get_file(),'r') as hf:
            # Printing attributes
            print('HDF5 file attributes:')
            for item in hf.attrs.keys():
                print('\t' + item + ":", hf.attrs[item])

    def get_allowed_h5exts(self):
        return self.__h5exts


class H5PICFile(H5File):
    def __init__(self, file, h5ftype=None):
        H5File.__init__(self, file, 'r')

        # Grid 3D types in filenames:
        self.__g3dtypes = ['density', 'field', 'current']
        self.__g3dsubgrid_str = 'subgrid'
        # RAW types in filenames:
        self.__rawtypes = ['raw']
        self.__n_time_chars = 8

        self.__h5ftype = h5ftype

    # Returning boolean: if file name contains 'raw'
    def is_g3d_file(self, fname):
        return any((mq in fname) for mq in self.__g3dtypes)

    # Returning boolean: if file name contains 'raw'
    def is_raw_file(self, fname):
        return any((mq in fname) for mq in self.__rawtypes)

    # Returning boolean:  if file extension is hdf5 extension and
    #                     if file name contains name of grid quantity
    def is_h5g3d_file(self):
        fname, fext = os.path.splitext(self.get_file())
        return self.is_h5_file(fext=fext) and self.is_g3d_file(fname=fname)

    # Returning boolean:  if file extension is hdf5 extension and
    #                     if file name contains 'raw'
    def is_h5raw_file(self):
        fname, fext = os.path.splitext(self.get_file())
        return self.is_h5_file(fext=fext) and self.is_raw_file(fname=fname)

    def fcheck(self):
        if self.__h5ftype == 'g3d':
            return self.is_h5g3d_file()
        elif self.__h5ftype == 'raw':       
            return self.is_h5raw_file()
        else:
            print('Error: No file type specified ["g3d", "raw"]!')
            sys.exit(1)

    def get_filename_time(self):
        name_w_time = os.path.splitext(self.get_filename())[0]
        stridx = [m.start() for m in re.finditer('_', name_w_time)][-1]
        return float(name_w_time[(stridx+1):])

    def get_filename_wo_time(self):
        name_w_time = os.path.splitext(self.get_filename())[0]
        stridx = [m.start() for m in re.finditer('_', name_w_time)][-1]
        name_wo_time = name_w_time[0:stridx]
        return name_wo_time

    def is_subgrid(self):
        fname = os.path.splitext(self.get_filename())[0]
        if (self.__g3dsubgrid_str in fname):
            return True
        else:
            return False


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
                                'p3':'p3',
                                'ipart' : 'ipart',
                                'iproc' : 'iproc'
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



class HiFile(H5Keys, H5PICFile):
    def __init__(self, file):
        H5Keys.__init__(self, 'hipace')
        H5PICFile.__init__(self, file)
        self.read_attrs()

    def read_attrs(self):
        # Reading attributes
        with h5py.File(self.get_file(),'r') as hf:
            self.__nx = hf.attrs[   self.get_attrkey('nx') ]
            self.__xmin = hf.attrs[ self.get_attrkey('xmin') ]
            self.__xmax = hf.attrs[ self.get_attrkey('xmax') ]
            self.__time = hf.attrs[ self.get_attrkey('time') ]
            self.__dt = hf.attrs[   self.get_attrkey('dt') ]
            type_bytes = hf.attrs[ self.get_attrkey('type') ]
            name_bytes = hf.attrs[ self.get_attrkey('name') ]

            if isinstance(type_bytes, str):
                self.__type = type_bytes
            else:
                self.__type = type_bytes[0].decode('UTF-8')
            if isinstance(name_bytes, str):
                self.__name = name_bytes
            else:
                self.__name = name_bytes[0].decode('UTF-8')

    def get_x_arr(self,dim):
        return np.linspace(self.__xmin[dim],self.__xmax[dim],self.__nx[dim])

    def get_z_arr(self):
        return np.linspace(self.__time+self.__xmin[0],self.__time+self.__xmax[0],self.__nx[0])

    def get_zeta_arr(self):
        return np.linspace(self.__xmin[0],self.__xmax[0],self.__nx[0])

    def get_xi_arr(self):
        return np.linspace(-self.__xmin[0],-self.__xmax[0],self.__nx[0])

    def get_nt(self):
        return round(self.__time/self.__dt)

    def get_time(self):
        return self.__time

    def get_nx(self,dim):
        return self.__nx[dim]        

    def get_xmin(self,dim):
        return self.__xmin[dim]

    def get_xmax(self,dim):
        return self.__xmax[dim]

    def get_dx(self,dim):
        return (self.__xmax[dim]-self.__xmin[dim])/self.__nx[dim]

    def get_dt(self):
        return self.__dt

    def get_name(self):
        return self.__name

    def get_type(self):
        return self.__type



class HiRAW(HiFile):
    def __init__(self, file):
        HiFile.__init__(self, file)
        self.__npart = 0
        self.__data_is_read = False

        self.__coord2idx = {
            "x1": 0,
            "x2": 1,
            "x3": 2,
            "p1": 3,
            "p2": 4,
            "p3": 5,
            "q": 6,
            "iproc": 7,
            "ipart": 8 }

    def __save_coord(self,coord, nparray):
        self.__part[self.__coord2idx[coord],:] = nparray

    def get(self,coord=''):
        if self.__data_is_read:
            if coord=='':
                return self.__part
            else:
                return self.__part[self.__coord2idx[coord]]
        else:
            print('Error: Data not yet read!')
            sys.exit(1) 

    def read_data(self, verbose=True):
        if verbose:
            print('Reading data: "%s" ...' % self.get_file())        
        with h5py.File(self.get_file(),'r') as hf:
            
            ncoords = 7

            # Get number of particles
            dsetq = hf[self.get_rawkey('q')]
            npart = (dsetq.shape)[0]

            if self.get_rawkey('ipart') in hf.keys() \
                and self.get_rawkey('iproc') in hf.keys():
                ncoords += 2

            # Allocating array
            self.__part = np.zeros((ncoords, npart), dtype=np.float32)

            # Reading datasets
            keys = list(self.__coord2idx.keys())
            for i in range(0,ncoords):
                self.__save_coord(  keys[i],\
                                    np.array(hf.get(  self.get_rawkey(keys[i]) )))

            self.__npart = npart

        self.__data_is_read = True
        
        if verbose:
            print('Data is read.') 

    def select_by_idx(self,idx):
        self.__part = self.__part[:,idx]

    def get_npart(self):
        self.__npart = (self.__part.shape)[1]
        return self.__npart

    def select_zeta_range(self, zeta_range, verbose=True):
        if not self.__data_is_read:
            self.read_data()
        if zeta_range != [] and len(zeta_range) == 2:
            idx = np.nonzero((self.get('x1') >= zeta_range[0]) & (self.get('x1') < zeta_range[1]))[0]
            self.__npart = np.size(idx)
            if verbose:
                print('%i particles in selected range [%0.2f, %0.2f]' \
                    % (self.__npart,zeta_range[0],zeta_range[1]))
            self.__part = self.__part[:,idx]


#### 3D-grid
class Grid3d(HiFile):
    def __init__(self, file):
        HiFile.__init__(self, file)
        self.is_read_3D = False
        self.read_attrs()
        with h5py.File(self.get_file(),'r') as hf:
            self.h5keys = list(hf.keys())
            if len(self.h5keys) == 1:
                self.dsetkey = self.h5keys[0]
            elif len(self.h5keys) == 0:
                print('Error:\tHDF5 file "%s" does not '
                  'contain any dataset!' %(self.get_file()) )
                sys.exit()
            elif len(self.h5keys) > 1:
                print('Error:\tHDF5 file "%s" contains more '
                  'than one dataset!' %(self.get_file()) )
                sys.exit()

        # self.fid = h5py.h5f.open(self.get_file().encode())
        # self.dset = h5py.h5d.open(self.fid, self.dsetkey.encode())

    def _read_3D(self):
        # read 3D, making sure data is read only once
        if self.is_read_3D == False:
            with h5py.File(self.get_file(),'r') as hf:
                # Reading dataset (here not caring how dataset is called)
                self.data3d = hf[self.dsetkey][()]
                self.is_read_3D == True
        return self.data3d

    def _read_2D(self, i0=None, i1=None, i2=None, navg=None):
        with h5py.File(self.get_file(),'r') as hf:
        # Reading dataset (here not caring how dataset is called)
            if i0!=None and i1==None and i2==None:
                data2d = hf[self.dsetkey][i0,:,:]
                if navg != None:
                    navgrg = int(navg/2)
                    norm = 0
                    for j0 in range(-navgrg,navgrg+1):
                        data2d += hf[self.dsetkey][i0+j0,:,:]
                        norm += 1
                    data2d = (data2d - hf[self.dsetkey][i0,:,:])/norm                 
            elif i0==None and i1!=None and i2==None:
                data2d = hf[self.dsetkey][:,i1,:]
                if navg != None:
                    navgrg = int(navg/2)
                    norm = 0
                    for j1 in range(-navgrg,navgrg+1):
                        data2d += hf[self.dsetkey][:,i1+j1,:]
                        norm += 1
                    data2d = (data2d - hf[self.dsetkey][:,i1,:])/norm                   
            elif i0==None and i1==None and i2!=None:
                data2d = hf[self.dsetkey][:,:,i2]
                if navg != None:
                    navgrg = int(navg/2)
                    norm = 0
                    for j2 in range(-navgrg,navgrg+1):
                        data2d += hf[self.dsetkey][:,:,i2+j2]
                        norm += 1
                    data2d = (data2d - hf[self.dsetkey][:,:,i2])/norm                 
            else:
                print('Error:\tExactly one index must '
                  'be provided for HDF slice read in!')
                sys.exit()
        return(data2d)

    def _read_1D(self, i0=None, i1=None, i2=None, navg=None):
        with h5py.File(self.get_file(),'r') as hf:
            # Reading dataset (here not caring how dataset is called)
            if i0!=None and i1!=None and i2==None:
                data1d = hf[self.dsetkey][i0,i1,:]
                if navg != None:
                    navgrg = int(navg/2)
                    norm = 0
                    for j0 in range(-navgrg,navgrg+1):
                        for j1 in range(-navgrg,navgrg+1):
                            data1d += hf[self.dsetkey][i0+j0,i1+j1,:]
                            norm += 1
                    data1d = (data1d - hf[self.dsetkey][i0,i1,:])/norm                 
            elif i0!=None and i1==None and i2!=None:
                data1d = hf[self.dsetkey][i0,:,i2]
                if navg != None:
                    navgrg = int(navg/2)
                    norm = 0
                    for j0 in range(-navgrg,navgrg+1):
                        for j2 in range(-navgrg,navgrg+1):
                            data1d += hf[self.dsetkey][i0+j0,:,i2+j2]
                            norm += 1
                    data1d = (data1d - hf[self.dsetkey][i0,:,i2])/norm             
            elif i0==None and i1!=None and i2!=None:
                data1d = hf[self.dsetkey][:,i1,i2]
                if navg != None:
                    navgrg = int(navg/2)
                    norm = 0
                    for j1 in range(-navgrg,navgrg+1):
                        for j2 in range(-navgrg,navgrg+1):
                            data1d += hf[self.dsetkey][:,i1+j1,i2+j2]
                            norm += 1
                    data1d = (data1d - hf[self.dsetkey][:,i1,i2])/norm                
            else:
                print('Error:\tExactly two indices must '
                  'be provided for HDF line read in!')
                sys.exit()
        return(data1d)

    def _read_0D(self, i0=None, i1=None, i2=None, navg=None):
        with h5py.File(self.get_file(),'r') as hf:
            # Reading dataset (here not caring how dataset is called)
            if i0!=None and i1!=None and i2!=None:
                data0d = hf[self.dsetkey][i0,i1,i2]
                if navg != None:
                    navgrg = int(navg/2)
                    norm = 0
                    for j0 in range(-navgrg,navgrg+1):
                        for j1 in range(-navgrg,navgrg+1):
                            for j2 in range(-navgrg,navgrg+1):
                                data0d += hf[self.dsetkey][i0+j0,i1+j1,i2+j2]
                                norm += 1
                    data0d = (data0d - hf[self.dsetkey][i0,i1,i2])/norm

            else:
                print('Error:\tExactly three indices must '
                  'be provided for HDF point read in!')
                sys.exit(1)
        return(data0d)

    def __n_none(self, arg0, arg1, arg2):
        return sum([arg0 == None, arg1 == None, arg2 == None])

    def read(self, 
             i0=None, 
             i1=None, 
             i2=None,
             x0=None,
             x1=None,
             x2=None,
             gradax=None,
             navg=None):
        """Read 3D data and return 0D, 1D, 2D or 3D gradient array.
        Args:
            i0 (int, optional): Index along axis0.
            i1 (int, optional): Index along axis1.
            i2 (int, optional): Index along axis2.
            x0 (float, optional): Position along axis0.
            x1 (float, optional): Position along axis1.
            x2 (float, optional): Position along axis2.
            ax (int, optional): Axis along which gradient is formed.
            navg (int, optional): Numer of gridpoints to be averaged for.                    

        Returns:
            ndarray: 0D, 1D, 2D or 3D numpy array
        """        
        if gradax != None:
            return self.read_grad(gradax,i0=i0,i1=i1,i2=i2,x0=x0,x1=x1,x2=x2)
        else:
            if self.__n_none(x0, x1, x2) == 3:
                if self.__n_none(i0, i1, i2) == 3:
                    return self._read_3D()
                elif self.__n_none(i0, i1, i2) == 2:
                    return self._read_2D(i0=i0, i1=i1, i2=i2, navg=navg) 
                elif self.__n_none(i0, i1, i2) == 1:
                    return self._read_1D(i0=i0, i1=i1, i2=i2, navg=navg)             
                elif self.__n_none(i0, i1, i2) == 0:
                    return self._read_0D(i0=i0, i1=i1, i2=i2, navg=navg)
            elif self.__n_none(i0, i1, i2) == 3:
                if x0 != None:
                    i0 = np.argmin(np.abs(self.get_x_arr(0)-x0))
                if x1 != None:
                    i1 = np.argmin(np.abs(self.get_x_arr(1)-x1))
                if x2 != None:
                    i2 = np.argmin(np.abs(self.get_x_arr(2)-x2))
                return self.read(i0=i0, i1=i1, i2=i2, navg=navg)                    
            else:
                print('Error:\tMixed indices and positions '
                  'not allowed for grid 3d read in!')
                sys.exit(1)                                  

    def read_integrate(self, 
                       ax0=False, 
                       ax1=False, 
                       ax2=False):
        # read and integrate along specified axes 
        data = self._read_3D()
        axtup = ()
        dx = 1
        if ax0:
            axtup += (0,)
            dx *= self.get_dx(0)
        if ax1:
            axtup += (1,)
            dx *= self.get_dx(1)
        if ax2:
            axtup += (2,)
            dx *= self.get_dx(2)           

        return np.sum(data,axis=axtup) * dx

    def __ax_in_plane(self, ax, i0, i1, i2):
        return ((ax+1) in (np.multiply([i0 != None, i1 != None, i2 != None],[1, 2, 3])))


    def read_grad(  self, 
                    ax,
                    i0=None, 
                    i1=None, 
                    i2=None,
                    x0=None,
                    x1=None,
                    x2=None):

        """Read 3D data, form derivative along given axis 
        and return 0D, 1D, 2D or 3D gradient array.
        
        Args:
            ax (int): Axis along which gradient is formed.

            dim (int): Dimension along which gradient is formed.
            i0 (int, optional): Index along axis0
            i1 (int, optional): Index along axis1
            i2 (int, optional): Index along axis2
            x0 (float, optional): Position along axis0
            x1 (float, optional): Position along axis1
            x2 (float, optional): Position along axis2

        Returns:
            ndarray: 0D, 1D, 2D or 3D gradient array
        """

        grad = np.gradient(self._read_3D(),axis=ax)/self.get_dx(ax)

        # read and differentiate along specified axes 
        if self.__n_none(x0, x1, x2) == 3:
            if self.__n_none(i0, i1, i2) == 3:
                out = grad
            elif self.__n_none(i0, i1, i2) == 2:
                if i0 != None:
                    out = grad[i0,:,:]
                elif i1 != None:
                    out = grad[:,i1,:]
                else:
                    out = grad[:,:,i2]
            elif self.__n_none(i0, i1, i2) == 1:
                if i0 == None:
                    out = grad[:,i1,i2]
                elif i1 == None:
                    out = grad[i0,:,i2]
                else:
                    out = grad[i0,i1,:]            
            elif self.__n_none(i0, i1, i2) == 0:
                out = grad[i0,i1,i2]
        elif self.__n_none(i0, i1, i2) == 3:
            if x0 != None:
                i0 = np.argmin(np.abs(self.get_x_arr(0)-x0))
            if x1 != None:
                i1 = np.argmin(np.abs(self.get_x_arr(1)-x1))
            if x2 != None:
                i2 = np.argmin(np.abs(self.get_x_arr(2)-x2))
            out = self.read_grad(ax,i0=i0, i1=i1, i2=i2)                    
        else:
            print('Error:\tMixed indices and positions '
              'not allowed for grid 3d read in!')
            sys.exit(1)

        return out

    def read_avgx(self, dim, ax0=False, 
                             ax1=False, 
                             ax2=False):
        # read and find avg position in dim and integrate along 
        # specified axes 
        data = self._read_3D()            

        axtup = ()
        dx = 1
        if dim == 0:
            sumx = np.dot( np.transpose(data,(1,2,0)),\
                           self.get_x_arr(0) )
            if ax0:
                print('Error:\tIntegration dimension must be '
                        'different from averaging dimension!')
                sys.exit(1)  
            if ax1:
                axtup += (0,)
                dx *= self.get_dx(1)
            if ax2:
                axtup += (1,)
                dx *= self.get_dx(2)  
                                            
        elif dim == 1:
            sumx = np.dot( np.transpose(data,(0,2,1)),\
                           self.get_x_arr(1) ) 
            if ax0:
                axtup += (0,)
                dx *= self.get_dx(0)
            if ax1:
                print('Error:\tIntegration dimension must be '
                        'different from averaging dimension!')
                sys.exit(1)                  
            if ax2:
                axtup += (1,)
                dx *= self.get_dx(2) 

        elif dim == 2:
            sumx = np.dot( np.transpose(data,(0,1,2)),\
                           self.get_x_arr(2) )
            if ax0:
                axtup += (0,)
                dx *= self.get_dx(0)
            if ax1:
                axtup += (1,)
                dx *= self.get_dx(1)                 
            if ax2:
                print('Error:\tIntegration dimension must be '
                        'different from averaging dimension!')
                sys.exit(1)                              
        else:
            print('Error:\tdim must be 0, 1 or 2!')
            sys.exit(1)                  
        
        norm = np.sum(data,axis=dim)
        norm[np.where(norm == 0)] = 1.0

        return np.divide(np.sum(sumx,axis=axtup),np.sum(norm,axis=axtup))

#### 3D-grid
class Grid2d(HiFile):
    def __init__(self, file):
        HiFile.__init__(self, file)
        self.read_attrs()
        with h5py.File(self.get_file(),'r') as hf:
            self.h5keys = list(hf.keys())
            if len(self.h5keys) == 1:
                self.dsetkey = self.h5keys[0]
            elif len(self.h5keys) == 0:
                print('Error:\tHDF5 file "%s" does not '
                  'contain any dataset!' %(self.get_file() ) )
                sys.exit()
            elif len(self.h5keys) > 1:
                print('Error:\tHDF5 file "%s" contains more '
                  'than one dataset!' %(self.get_file() ) )
                sys.exit()

        # self.fid = h5py.h5f.open(self.file.encode())
        # self.dset = h5py.h5d.open(self.fid, self.dsetkey.encode())

    def _read_3D(self):
        print('Error:\tGrid2D can not read in 3D dataset!' )
        sys.exit()
        return 


    def _read_2D(self, i0=None, i1=None, i2=None):
        with h5py.File(self.get_file(),'r') as hf:
        # Reading dataset (here not caring how dataset is called)
            
            data2d = hf[self.dsetkey][()]
            # elif i0==None and i1!=None and i2==None:
            #     data2d = hf[self.dsetkey][:,i1]
            # elif i0==None and i1==None and i2!=None:
            #     data2d = hf[self.dsetkey][:,i2]
        return(data2d)

    def _read_1D(self, i0=None, i1=None, i2=None):
        with h5py.File(self.get_file(),'r') as hf:
            # Reading dataset (here not caring how dataset is called)
            if i0!=None and i1==None:
                data1d = hf[self.dsetkey][i0,:]
            elif i0==None and i1!=None:
                data1d = hf[self.dsetkey][:,i1]
            elif i0==None and i2!=None:
                data1d = hf[self.dsetkey][:,i2]
                print('Warning: Wrong index (only 2d array available). \n'
                        'Switched to x1 instead of x2!')
            else:
                print('Error:\t Something went wrong with the reading of the line out')
                sys.exit()
        return(data1d)

    def _read_0D(self, i0=None, i1=None, i2=None):
        with h5py.File(self.get_file(),'r') as hf:
            # Reading dataset (here not caring how dataset is called)
            if i0!=None and i1!=None and i2!=None:
                data0d = hf[self.dsetkey][i0,i1,i2]
            else:
                print('Error:\tExactly three indices must '
                  'be provided for HDF point read in!')
                sys.exit(1)
        return(data0d)

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
                return self._read_2D()
            elif self.__n_none(i0, i1, i2) == 2:
                return self._read_1D(i0=i0, i1=i1, i2=i2) 
            elif self.__n_none(i0, i1, i2) == 1:
                return self._read_0D(i0=i0, i1=i1, i2=i2)             
            elif self.__n_none(i0, i1, i2) == 0:
                print('Error:\tToo many indices for 2D array! ')
                sys.exit(1)    
                return self._read_0D(i0=i0, i1=i1, i2=i2)
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

    def read_integrate(self, 
                       ax0=False, 
                       ax1=False, 
                       ax2=False):
        # read and integrate along specified axes 
        data = self._read_3D()
        axtup = ()
        dx = 1
        if ax0:
            axtup += (0,)
            dx *= self.get_dx(0)
        if ax1:
            axtup += (1,)
            dx *= self.get_dx(1)
        if ax2:
            axtup += (2,)
            dx *= self.get_dx(2)           

        return np.sum(data,axis=axtup) * dx


    def read_avgx(self, dim, ax0=False, 
                             ax1=False, 
                             ax2=False):
        # read and find avg position in dim and integrate along 
        # specified axes 
        data = self._read_3D()            

        axtup = ()
        dx = 1
        if dim == 0:
            sumx = np.dot( np.transpose(data,(1,2,0)),\
                           self.get_x_arr(0) )
            if ax0:
                print('Error:\tIntegration dimension must be '
                        'different from averaging dimension!')
                sys.exit(1)  
            if ax1:
                axtup += (0,)
                dx *= self.get_dx(1)
            if ax2:
                axtup += (1,)
                dx *= self.get_dx(2)  
                                            
        elif dim == 1:
            sumx = np.dot( np.transpose(data,(0,2,1)),\
                           self.get_x_arr(1) ) 
            if ax0:
                axtup += (0,)
                dx *= self.get_dx(0)
            if ax1:
                print('Error:\tIntegration dimension must be '
                        'different from averaging dimension!')
                sys.exit(1)                  
            if ax2:
                axtup += (1,)
                dx *= self.get_dx(2) 

        elif dim == 2:
            sumx = np.dot( np.transpose(data,(0,1,2)),\
                           self.get_x_arr(2) )
            if ax0:
                axtup += (0,)
                dx *= self.get_dx(0)
            if ax1:
                axtup += (1,)
                dx *= self.get_dx(1)                 
            if ax2:
                print('Error:\tIntegration dimension must be '
                        'different from averaging dimension!')
                sys.exit(1)                              
        else:
            print('Error:\tdim must be 0, 1 or 2!')
            sys.exit(1)                  
        
        norm = np.sum(data,axis=dim)
        norm[np.where(norm == 0)] = 1.0

        return np.divide(np.sum(sumx,axis=axtup),np.sum(norm,axis=axtup))

class SliceMoms(H5File):
    def __init__(self):
        self.__time_array_name = 'time_array'
        self.__zeta_array_name = 'zeta_array'

        self.__mom_order_h5key = 'mom_order'
        self.__xterms_h5key = 'xterms'

        self.__order = 2
        self.__Nzeta = None
        self.__Nt = None
        self.__arrays_allocated = False
        self.__data_is_read = False
        self.__with_2nd_order_xterms = False

        self.__set_mom2idx_dict()

    def __set_mom2idx_dict(self):
        self.__mom2idx_dict = {}

        self.__mom2idx_dict['0'] = {}
        self.__mom2idx_dict['0']['000000'] = 0
        
        self.__mom2idx_dict['1'] = {}
        self.__mom2idx_dict['1']['100000'] = 0
        self.__mom2idx_dict['1']['010000'] = 1
        self.__mom2idx_dict['1']['001000'] = 2
        self.__mom2idx_dict['1']['000100'] = 3
        self.__mom2idx_dict['1']['000010'] = 4
        self.__mom2idx_dict['1']['000001'] = 5

        self.__mom2idx_dict['2'] = {}
        self.__mom2idx_dict['2']['200000'] = 0
        self.__mom2idx_dict['2']['020000'] = 1
        self.__mom2idx_dict['2']['002000'] = 2
        self.__mom2idx_dict['2']['000200'] = 3
        self.__mom2idx_dict['2']['000020'] = 4
        self.__mom2idx_dict['2']['000002'] = 5
        self.__mom2idx_dict['2']['100100'] = 6
        self.__mom2idx_dict['2']['010010'] = 7
        self.__mom2idx_dict['2']['001001'] = 8

        self.__mom2idx_dict['2x'] = {}
        self.__mom2idx_dict['2x']['200000'] = 0
        self.__mom2idx_dict['2x']['020000'] = 1
        self.__mom2idx_dict['2x']['002000'] = 2
        self.__mom2idx_dict['2x']['000200'] = 3
        self.__mom2idx_dict['2x']['000020'] = 4
        self.__mom2idx_dict['2x']['000002'] = 5
        self.__mom2idx_dict['2x']['110000'] = 6
        self.__mom2idx_dict['2x']['101000'] = 7
        self.__mom2idx_dict['2x']['100100'] = 8
        self.__mom2idx_dict['2x']['100010'] = 9
        self.__mom2idx_dict['2x']['100001'] = 10
        self.__mom2idx_dict['2x']['011000'] = 11
        self.__mom2idx_dict['2x']['010100'] = 12
        self.__mom2idx_dict['2x']['010010'] = 13
        self.__mom2idx_dict['2x']['010001'] = 14
        self.__mom2idx_dict['2x']['001100'] = 15
        self.__mom2idx_dict['2x']['001010'] = 16
        self.__mom2idx_dict['2x']['001001'] = 17
        self.__mom2idx_dict['2x']['000110'] = 18
        self.__mom2idx_dict['2x']['000101'] = 19
        self.__mom2idx_dict['2x']['000011'] = 20

        self.__mom2idx_dict['3'] = {}
        self.__mom2idx_dict['3']['300000'] = 0
        self.__mom2idx_dict['3']['030000'] = 1
        self.__mom2idx_dict['3']['003000'] = 2
        self.__mom2idx_dict['3']['000300'] = 3
        self.__mom2idx_dict['3']['000030'] = 4
        self.__mom2idx_dict['3']['000003'] = 5
        self.__mom2idx_dict['3']['200100'] = 6
        self.__mom2idx_dict['3']['100200'] = 7
        self.__mom2idx_dict['3']['020010'] = 8
        self.__mom2idx_dict['3']['010020'] = 9
        self.__mom2idx_dict['3']['002001'] = 10
        self.__mom2idx_dict['3']['001002'] = 11

        self.__mom2idx_dict['4'] = {}
        self.__mom2idx_dict['4']['400000'] = 0
        self.__mom2idx_dict['4']['040000'] = 1
        self.__mom2idx_dict['4']['004000'] = 2
        self.__mom2idx_dict['4']['000400'] = 3
        self.__mom2idx_dict['4']['000040'] = 4
        self.__mom2idx_dict['4']['000004'] = 5
        self.__mom2idx_dict['4']['200200'] = 6
        self.__mom2idx_dict['4']['020020'] = 7
        self.__mom2idx_dict['4']['002002'] = 8   

    def get_h5order(self, file):
        with h5py.File(file,'r') as hf:
            if self.__mom_order_h5key in hf.attrs.keys():
                h5order  = hf.attrs[self.__mom_order_h5key]
            else:
                h5order = 2
        return h5order

    def get_if_xterms(self):
        return self.__with_2nd_order_xterms

    def get_if_h5xterms(self, file):
        with h5py.File(file,'r') as hf:
            if self.__xterms_h5key in hf.attrs.keys():
                xterms  = hf.attrs[self.__xterms_h5key]
            else:
                xterms = False

        if xterms:
            print('With x-terms!')            

        return xterms

    def __get_n_moments(self, order, with_2nd_order_xterms = False):
        n_moments = 0
        for i in range(order+1):
            if order == 2 and with_2nd_order_xterms == True:
                key = '2x'
            else:
                key = '%d' % i
            n_moments += len(self.__mom2idx_dict[key])
        return n_moments

    def alloc(self, Nzeta, Nt, order = 2, with_2nd_order_xterms = False):
        self.__order = order
        self.__with_2nd_order_xterms = with_2nd_order_xterms
        self.__Nzeta = Nzeta
        self.__Nt = Nt

        self.__zeta_array = np.zeros(Nt, dtype=np.float32)   
        self.__time_array = np.zeros(Nt, dtype=np.float32)

        if with_2nd_order_xterms:
            print('With x-terms!')          

        n_moments = self.__get_n_moments(order, with_2nd_order_xterms)

        print('Number of moments: %d' % n_moments)

        if order <= 4:
            self.__mom = np.zeros((n_moments, Nt, Nzeta), dtype=np.float32)
        else:
            print('Error:\tOnly moment orders up to 4 implemented!')
            sys.exit(1)              

        self.__arrays_allocated = True

    def read(self, file, order = None, verbose = True):
        H5File.__init__(self, file, 'r')
        if not self.is_h5_file():
            print('Error:\tFile is not an HDF5 file!')
            sys.exit(1)      

        if verbose:
            print('Reading %s ...' % file)

        h5order = self.get_h5order(file)
        self.__with_2nd_order_xterms = self.get_if_h5xterms(file)

        with h5py.File(self.get_file(),'r') as hf:
            if order == None:
                order = h5order
            else:
                if order > h5order:
                    print('Error:\tSpecified order is greater than moment order in hdf5 file!')
                    sys.exit(1)     
            self.__order = order

            # Reading datasets
            self.__time_array = np.array(hf.get( self.__time_array_name ))
            self.__zeta_array = np.array(hf.get( self.__zeta_array_name ))

            if verbose:
                print('Reading moments. Order: %i' % self.__order)
            self.__mom = np.array(hf.get( 'mom' ))

        self.__data_is_read = True

    def get_zeta_array(self):
        if not self.__data_is_read:
            print('Error:\tArrays have not been read!')
            sys.exit(1)   
        return self.__zeta_array

    def get_time_array(self):
        if not self.__data_is_read:
            print('Error:\tArrays have not been read!')
            sys.exit(1)   
        return self.__time_array


    def __orders_to_idx(self, x1, x2, x3, p1, p2, p3, err_ignore=False):  

        mom_str = '%d%d%d%d%d%d' % (x1,x2,x3,p1,p2,p3)

        order = sum([x1,x2,x3,p1,p2,p3])

        offset = self.__get_n_moments(order-1, self.__with_2nd_order_xterms)

        if order == 2 and self.__with_2nd_order_xterms == True:
            key = '2x'
        else:
            key = '%d' % order        

        idx = offset + self.__mom2idx_dict[key][mom_str]

        return idx
  

    def get(self, x=0, y=0, z=0, px=0, py=0, pz=0, orders = None):
        if not self.__data_is_read:
            print('Error:\tArrays have not been read!')
            sys.exit(1)

        if orders != None:
            if len(orders) == 6:
                # Overwriting single orders
                x = orders[0]
                y = orders[1]
                z = orders[2]
                px = orders[3]
                py = orders[4]
                pz = orders[5]
            else:
                print('Error:\tProvided coordinates must have length 6!')
                sys.exit(1)                

        idx = self.__orders_to_idx( x1=z, x2=x, x3=y, \
                                    p1=pz, p2=px, p3=py)

        return self.__mom[idx,:,:]

    def get_charge(self):
        return self.get()       

    def get_noncentral(self, x=0, y=0, z=0, px=0, py=0, pz=0):

        centralization = np.ones(self.get().shape)

        orders = [x,y,z,px,py,pz]
        if sum(orders) == 2:    
            for i in range(len(orders)):
                uniorders = [0,0,0,0,0,0]
                uniorders[i] = 1
                centralization = np.multiply(centralization,\
                    np.power(self.get(orders = uniorders),orders[i]))

            noncentral = self.get(x=x, y=y, z=z, px=px, py=py, pz=pz) + centralization
                    
        else:
            print('Error:\tCalculation of noncentral '
                    'moments with order > 2 not yet implemented!')
            sys.exit(1)            

        return noncentral


    def project(self,quantity):
        beam_charge = np.sum(self.get_charge(), axis=1)
        return np.divide(
                np.sum( \
                    np.multiply(quantity,
                                self.get_charge()),\
                    axis=1),\
                beam_charge)        


    def get_proj(self, x=0, y=0, z=0, px=0, py=0, pz=0):

        orders = [x,y,z,px,py,pz]

        if sum(orders) == 1:

            # Compute simple projection of first order moment
            proj = self.project(
                    self.get(   x=x, \
                                y=y, \
                                z=z, \
                                px=px, \
                                py=py, \
                                pz=pz))
        elif sum(orders) == 2:
            
            noncentral_proj = self.project(
                self.get_noncentral(x=x, \
                                    y=y, \
                                    z=z, \
                                    px=px, \
                                    py=py, \
                                    pz=pz) )

            proj_centralization = np.ones(noncentral_proj.shape)

            for i in range(len(orders)):
                uniorders = [0,0,0,0,0,0]
                uniorders[i] = 1
                proj_centralization = \
                    np.multiply(proj_centralization,\
                        np.power(\
                                self.project(\
                                    self.get(orders = uniorders)),\
                                orders[i]))

            # Compute 2nd order central moment
            # according to <(x-<x>)^2> = <x^2> - <x>^2
            # where <*> is the non-central average
            central_proj = noncentral_proj - proj_centralization

            if np.count_nonzero(orders) == 1:
                # Use abs to make sure floating point errors don't 
                #  lead to negative central squared moments
                proj = np.abs(central_proj)
            else:
                proj = central_proj
                    
        else:
            print('Error:\tCalculation of noncentral '
                    'moments with order > 2 not yet implemented!')
            sys.exit(1)  

        return proj

    def set_time(self,time,nt):
        if not self.__arrays_allocated:
            print('Error:\tArrays need to be allocated!')
            sys.exit(1)
        if nt >= self.__Nt:
            print('Error:\tnt exceeds Nt!')
            sys.exit(1)                    
        self.__time_array[nt] = time

    def set_zeta_array(self,nparray):
        if not self.__arrays_allocated:
            print('Error:\tArrays need to be allocated!')
            sys.exit(1)

        if np.shape(nparray) != (self.__Nzeta,):
            print('Error:\tArray does not match expected size (%i,)!' % self.__Nzeta)
            sys.exit(1)
        self.__zeta_array = nparray

    def set_at_nt(self, nparray, nt, x1=0, x2=0, x3=0, p1=0, p2=0, p3=0):
        if not self.__arrays_allocated:
            print('Error:\tArrays need to be allocated!')
            sys.exit(1)

        if type(nparray) is not np.ndarray:
            print('Error:\targument "nparray" must be numpy array!')
            sys.exit(1)

        if np.shape(nparray) != (self.__Nzeta,):
            print('Error:\tArray does not match expected size (%i,)!' % self.__Nzeta)
            sys.exit(1)

        if x1+x2+x3+p1+p2+p3 > self.__order:
            print('Error:\tAttempting to set moment with order ' \
                'higher than expected (%i vs. %i)!' % (x1+x2+x3+p1+p2+p3,self.__order))
            sys.exit(1)

        idx = self.__orders_to_idx(x1,x2,x3,p1,p2,p3)

        self.__mom[idx,nt,:] = nparray       

    def write(self, file, order = None):
        H5File.__init__(self, file, 'w')
        if order == None:
            order = self.__order

        if not self.is_h5_file():
            print('Error:\tFile name does not have an hdf5 ending!')
            sys.exit(1) 

        h5f = h5py.File(file, "w")

        h5f.attrs[self.__mom_order_h5key] = order
        h5f.attrs[self.__xterms_h5key] = self.__with_2nd_order_xterms

        if not self.__arrays_allocated:
            print('Error:\tArrays do not exist!')
            sys.exit(1)

        dset_zeta_array = h5f.create_dataset( self.__zeta_array_name, data = self.__zeta_array )
        dset_time_array = h5f.create_dataset( self.__time_array_name, data = self.__time_array )
        dset_mom = h5f.create_dataset(  "mom", data = self.__mom )

        h5f.close()

    def truncate_zeta_region(self, zeta_min, zeta_max, order = None, crossterms = False):

        if order == None:
            order = self.__order

        idx = np.nonzero(np.logical_and( self.__zeta_array >= zeta_min, \
            self.__zeta_array <= zeta_max ))[0]

        self.__zeta_array = self.__zeta_array[idx]
        self.__mom = self.__mom[:,:,idx]

    def shift_time(self, time_offset):
        if time_offset != 0.0:
            print('Shifting time array by: %f' % time_offset)
        self.__time_array += time_offset 


class H5Plot:
    def __init__(self):

        # Allowed hdf5 extensions:
        self.__h5exts = ['.h5','.hdf5']

        self.__xlab_key = u'xlab'
        self.__ylab_key = u'ylab'
        self.__xlab = u'$x$'
        self.__ylab = u'$y$'

        self.__lp_key = u'line_plots'
        self.__lp_label_key = u'label'  
        self.__lp_linestyle_key = u'linestyle'
        self.__lp_color_key = u'color'        
        self.__lp_x_key = u'x'
        self.__lp_y_key = u'y'

        self.__N_lp = 0
        self.__lp_X = []
        self.__lp_Y = []
        self.__lp_linestyles = []
        self.__lp_colors = []
        self.__lp_labels = []

    # Returning boolean: if file extension is hdf5 extension
    def is_h5_file(self, file):
        fname, fext = os.path.splitext(file)
        return any(fext == h5ext for h5ext in self.__h5exts)   

    def set_ax_labels(xlab, ylab):
        self.__xlab = np.string_(xlab)
        self.__ylab = np.string_(ylab)        

    def get_xlab(self):
        return self.__xlab 

    def get_ylab(self):
        return self.__ylab        

    def get_ax(self):
        return self.__ax 

    def inherit_matplotlib_line_plots(self, ax):
        self.__ax = ax
        self.__xlab = ax.xaxis.get_label_text()
        self.__ylab = ax.yaxis.get_label_text()
        self.__N_lp = len(ax.lines)
        for i in range(0,self.__N_lp):
            line = ax.lines[i]
            self.__lp_X.append(line.get_xdata())
            self.__lp_Y.append(line.get_ydata())
            self.__lp_labels.append(line.get_label())
            self.__lp_linestyles.append(line.get_linestyle())
            self.__lp_colors.append(line.get_color())

    def append_line_plot(self, x, y, label):
        self.__N_lp += 1
        self.__lp_X.append(x)
        self.__lp_Y.append(y)
        self.__lp_labels.append(np.string_(label))

    def write(self, path):
        h5f = h5py.File(path, "w")
        if self.__N_lp > 0:
            lp_grp = h5f.create_group(self.__lp_key)
            lp_grp.attrs[self.__xlab_key] = self.__xlab
            lp_grp.attrs[self.__ylab_key] = self.__ylab
            subgrps = []
            for i in range(0,self.__N_lp):
                subgrps.append(lp_grp.create_group(u'line_%02d' % i))
                subgrps[i].attrs[self.__lp_label_key] = self.__lp_labels[i]
                subgrps[i].attrs[self.__lp_linestyle_key] = self.__lp_linestyles[i]
                subgrps[i].attrs[self.__lp_color_key] = self.__lp_colors[i]
                subgrps[i].create_dataset( self.__lp_x_key, data = self.__lp_X[i])
                subgrps[i].create_dataset( self.__lp_y_key, data = self.__lp_Y[i])
        h5f.close()

    def read(self, path):
        self.__init__()
        if not self.is_h5_file(path):
            print('ERROR: Provided file %s is no hdf5 file!' % path)
            sys.exit(1)

        h5f = h5py.File(path, "r")
        self.__xlab = h5f[self.__lp_key].attrs[self.__xlab_key]
        self.__ylab = h5f[self.__lp_key].attrs[self.__ylab_key]
        subgrps = [] 
        for key in h5f[self.__lp_key].keys():
            self.__N_lp += 1
            i = self.__N_lp - 1
            subgrps.append(h5f[self.__lp_key].get(key))
            self.__lp_X.append(np.array(subgrps[i].get(self.__lp_x_key)))
            self.__lp_Y.append(np.array(subgrps[i].get(self.__lp_y_key)))
            self.__lp_labels.append(subgrps[i].attrs[self.__lp_label_key])
            self.__lp_linestyles.append(subgrps[i].attrs[self.__lp_linestyle_key])
            self.__lp_colors.append(subgrps[i].attrs[self.__lp_color_key])       
        h5f.close()

    def get_line_plots(self):
        return zip(self.__lp_X, self.__lp_Y, self.__lp_labels, self.__lp_linestyles, self.__lp_colors)

    def get_data(self, idx=0):
        return self.__lp_X[idx], self.__lp_Y[idx]      

class H5FList():
    def __init__(self, paths, h5ftype=None):
        self.__paths = paths
        self.__h5ftype = h5ftype
        self.__flist = None

    def get(self, verbose=True, stride=1):
        if not self.__paths:
            print('Error: No file provided!')
            sys.exit(1)

        if isinstance(self.__paths, list):
            # if 'paths' is a list of directories
            list_of_flists = []
            for path in self.__paths:
                list_of_flists.append(self.__get_individ(path, verbose))
            flist = [item for sublist in list_of_flists for item in sublist]
        elif isinstance(self.__paths, str):
            # if 'paths' is a single directory
            flist = self.__get_individ(self.__paths, verbose)

        # Alphabetically sorting list
        self.__flist = sorted(flist)
        return self.__flist[0::stride]

    def __get_individ(self, path, verbose):
        flist = []
        if os.path.isfile(path):
            file = path
            h5f = H5PICFile(file, h5ftype=self.__h5ftype)
            if h5f.fcheck():
                flist.append(file)
            else:
                if verbose: print('Skipping: ' + file)
        elif os.path.isdir(path):
            if verbose: print('"' + path + '"' + ' is a directory.')
            if verbose: 
                print('Processing all ' + self.__h5ftype + 
                      ' files in the provided directory.')
            for root, dirs, files in os.walk(path):
                for filename in files:
                    file = root + '/' + filename
                    h5f = H5PICFile(file, h5ftype=self.__h5ftype)
                    if h5f.fcheck():
                        flist.append(file)
                    else:
                       if verbose: print('Skipping: ' + file)
        elif not os.path.exists(path):
            print('Error: Provided path "%s" does not exist!' % path)
            sys.exit()
        else:
            print('Error: Provided path "%s" is neither a file nor a directory!' % path)
            sys.exit()
        return flist 

    def get_uniques(self, n_time_chars = 8):
        fnames = []
        if self.__flist == None:
            self.get()
        for f in self.__flist:
            h5f = H5PICFile(f)
            fnames.append(h5f.get_filename_wo_time())
        return list(set(fnames))

    def split_by_uniques(self, n_time_chars = 8):
        if self.__flist == None:
            self.get()
        uniques = self.get_uniques(n_time_chars=n_time_chars)

        # initialize and append to list of lists
        lofl = [[] for i in range(len(uniques))] 
        for i in range(len(uniques)):
            for f in self.__flist:
                if uniques[i] in os.path.split(f)[1]:
                    lofl[i].append(f)
        return lofl

    def get_paths(self):
        return self.__paths

    def get_h5ftype(self):
        return self.__h5ftype



def _localize_path_ifnabs(path):
    # If path is not an absolute path and not explicitly a local path,
    # assume path is local and prepend './'
    if path != None:
        if (not os.path.isabs(path)) and (not path[0] == '.'):
            path = './' + path
    else:
        path = None
    return path

def mkdirs_if_nexist(path, localize=True):
    """Function which recursively generates directories.

    Args:
        path (str): Path for which directories are generated.
    """
    
    if localize:
        path = _localize_path_ifnabs(path)

    return_path = path

    if path != None:
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

    return return_path
