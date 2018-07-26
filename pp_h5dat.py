#!/usr/bin/env python3
# pp_h5dat.py

import os
import sys
import numpy as np
import h5py


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
        self.__g3dsubgrid_str = 'subgrid'
        # RAW types in filenames:
        self.__rawtypes = ['raw']
        self.__n_time_chars = 6

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

    def get_filename(self):
        return os.path.split(self.file)[1]

    def get_filename_time(self):
        name_w_time = os.path.splitext(os.path.split(self.file)[1])[0]
        return float(name_w_time[-self.__n_time_chars:])

    def get_filename_wo_time(self):
        name_w_time = os.path.splitext(os.path.split(self.file)[1])[0]
        name_wo_time = name_w_time[0:-self.__n_time_chars]
        return name_wo_time

    def is_subgrid(self):
        fname = os.path.splitext(os.path.split(self.file)[1])[0]
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

    def get_time(self):
        return self.time

    def get_nx(self,dim):
        return self.nx[dim]        

    def get_xmin(self,dim):
        return self.xmin[dim]

    def get_xmax(self,dim):
        return self.xmax[dim]

    def get_dx(self,dim):
        return (self.xmax[dim]-self.xmin[dim])/self.nx[dim]

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
        self.__npart = 0
        self.__data_is_read = False

    def read_data(self, verbose=True):
        if verbose:
            print('Reading data...')        
        with h5py.File(self.file,'r') as hf:
            # Reading datasets
            self.x1 = np.array(hf.get(  self.get_rawkey('x1') ))
            self.x2 = np.array(hf.get(  self.get_rawkey('x2') ))
            self.x3 = np.array(hf.get(  self.get_rawkey('x3') ))
            self.q =  np.array(hf.get(  self.get_rawkey('q')  ))
            self.p1 = np.array(hf.get(  self.get_rawkey('p1') ))
            self.p2 = np.array(hf.get(  self.get_rawkey('p2') ))
            self.p3 = np.array(hf.get(  self.get_rawkey('p3') ))

            self.__npart = len(self.q)

        self.__data_is_read = True
        if verbose:
            print('Data is read.') 

    def get_npart(self):
        return self.__npart

    def select_zeta_range(self, zeta_range, verbose=True):
        if not self.__data_is_read:
            self.read_data()
        if zeta_range != [] and len(zeta_range) == 2:
            idx = np.nonzero((self.x1 >= zeta_range[0]) & (self.x1 < zeta_range[1]))
            self.__npart = np.size(idx)
            if verbose:
                print('%i particles in selected range [%0.2f, %0.2f]' % (self.__npart,zeta_range[0],zeta_range[1]))
            self.x1 = self.x1[idx]
            self.x2 = self.x2[idx]
            self.x3 = self.x3[idx]
            self.p1 = self.p1[idx]
            self.p2 = self.p2[idx]
            self.p3 = self.p3[idx]
            self.q = self.q[idx]

    # def get_var(self, var, idx = [], ifread = True):
    #     ret_array = []
    #     if not self.__data_is_read:
    #         if ifread:
    #             if var in self.__rawkeys:
    #                 with h5py.File(self.file,'r') as hf:
    #                     exec("self.%s = np.array(hf.get(  self.get_rawkey('%s') ))" % (var,var))
    #             else:
    #                 print('Error:\tVariable "%s" does not exist!' % var )
    #                 sys.exit()                     
    #         else:    
    #             print('Error:\tFile %s has not been read!' % (self.file) )
    #             sys.exit()
    #     exec('ret_array = self.%s' % var)  
    #     exec('print(self.%s)' % var)
    #     print(ret_array)
    #     if idx == []:
    #         return ret_array
    #     else:
    #         return ret_array[idx]

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
            data3d = hf[self.dsetkey][()]
        return data3d


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

    def read_integrate(self, 
                       ax0=False, 
                       ax1=False, 
                       ax2=False):
        # read and integrate along specified axes 
        data = self.read_3D()
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
        data = self.read_3D()            

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
        self.time_array_name = 'time_array'
        self.zeta_array_name = 'zeta_array'
        self.__order = -1

    def get_order(self):
        return self.__order 

    def alloc(self, Nzeta, Nt, order = 2):

        self.__order = order

        self.zeta_array = np.zeros(Nt, dtype=np.float32)   
        self.time_array = np.zeros(Nt, dtype=np.float32)   
        self.charge = np.zeros((Nt, Nzeta), dtype=np.float32)

        if order>0:
            self.avgx1 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx2 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx3 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgp1 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgp2 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgp3 = np.zeros((Nt, Nzeta), dtype=np.float32)

        if order>1:
            self.avgx1sq = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx2sq = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx3sq = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgp1sq = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgp2sq = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgp3sq = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx1p1 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx2p2 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx3p3 = np.zeros((Nt, Nzeta), dtype=np.float32)

            # crossterms:
            self.avgx1x2 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx3x1 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx2x3 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgp1p2 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgp3p1 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgp2p3 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx1p2 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx1p3 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx2p1 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx2p3 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx3p1 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx3p2 = np.zeros((Nt, Nzeta), dtype=np.float32)

        if order>2:
            self.avgx1cube = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx2cube = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx3cube = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgp1cube = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgp2cube = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgp3cube = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx1sqp1 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx2sqp2 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx3sqp3 = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx1p1sq = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx2p2sq = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx3p3sq = np.zeros((Nt, Nzeta), dtype=np.float32)

        if order>2:
            self.avgx1quar = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx2quar = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx3quar = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgp1quar = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgp2quar = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgp3quar = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx1sqp1sq = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx2sqp2sq = np.zeros((Nt, Nzeta), dtype=np.float32)
            self.avgx3sqp3sq = np.zeros((Nt, Nzeta), dtype=np.float32)

    def read(self, file, order = None, verbose = True):
        H5File.__init__(self, file)
        if not self.is_h5_file():
            print('Error:\tFile is not an HDF5 file!')
            sys.exit(1)      

        if verbose:
            print('Reading %s ...' % file)

        with h5py.File(self.file,'r') as hf:
            # TODO: read order of moments and if crossterms from attributes!!!

            if 'mom_order' in hf.keys():
                h5_mom_order = hf.attrs['mom_order']
            else:
                h5_mom_order = 2

            if order == None:
                order = h5_mom_order
            else:
                if order > h5_mom_order:
                    print('Error:\tSpecified order is greater than moment order in hdf5 file!')
                    sys.exit(1)     
            self.__order = order

            # Reading datasets
            self.time_array = np.array(hf.get( self.time_array_name ))
            self.zeta_array = np.array(hf.get( self.zeta_array_name ))
            self.charge = np.array(hf.get( 'charge' ))

            if verbose:
                print('Reading first order moments...')
            first_order = hf.get('first_order')
            self.avgx1 = np.array(first_order.get( 'avgx1' ))
            self.avgx2 = np.array(first_order.get( 'avgx2' ))
            self.avgx3 = np.array(first_order.get( 'avgx3' ))
            self.avgp1 = np.array(first_order.get( 'avgp1' ))
            self.avgp2 = np.array(first_order.get( 'avgp2' ))
            self.avgp3 = np.array(first_order.get( 'avgp3' ))

            if order > 1:
                if verbose:
                    print('Reading second order moments...')
                second_order = hf.get('second_order')
                self.avgx1sq = np.array(second_order.get( 'avgx1sq' ))
                self.avgx2sq = np.array(second_order.get( 'avgx2sq' ))
                self.avgx3sq = np.array(second_order.get( 'avgx3sq' ))
                self.avgp1sq = np.array(second_order.get( 'avgp1sq' ))
                self.avgp2sq = np.array(second_order.get( 'avgp2sq' ))
                self.avgp3sq = np.array(second_order.get( 'avgp3sq' ))
                self.avgx1p1 = np.array(second_order.get( 'avgx1p1' ))
                self.avgx2p2 = np.array(second_order.get( 'avgx2p2' ))
                self.avgx3p3 = np.array(second_order.get( 'avgx3p3' ))

            if order > 2:
                if verbose:
                    print('Reading third order moments...')
                third_order = hf.get('third_order')
                self.avgx1cube = np.array(third_order.get( 'avgx1cube' ))
                self.avgx2cube = np.array(third_order.get( 'avgx2cube' ))
                self.avgx3cube = np.array(third_order.get( 'avgx3cube' ))
                self.avgp1cube = np.array(third_order.get( 'avgp1cube' ))
                self.avgp2cube = np.array(third_order.get( 'avgp2cube' ))
                self.avgp3cube = np.array(third_order.get( 'avgp3cube' ))
                self.avgx1sqp1 = np.array(third_order.get( 'avgx1sqp1' ))
                self.avgx2sqp2 = np.array(third_order.get( 'avgx2sqp2' ))
                self.avgx3sqp3 = np.array(third_order.get( 'avgx3sqp3' ))
                self.avgx1p1sq = np.array(third_order.get( 'avgx1p1sq' ))
                self.avgx2p2sq = np.array(third_order.get( 'avgx2p2sq' ))
                self.avgx3p3sq = np.array(third_order.get( 'avgx3p3sq' ))                

                # crossterms:
                # ...

            if order > 3:
                if verbose:
                        print('Reading fourth order moments...')
                fourth_order = hf.get('fourth_order')
                self.avgx1quar = np.array(fourth_order.get( 'avgx1quar' ))
                self.avgx2quar = np.array(fourth_order.get( 'avgx2quar' ))
                self.avgx3quar = np.array(fourth_order.get( 'avgx3quar' ))
                self.avgp1quar = np.array(fourth_order.get( 'avgp1quar' ))
                self.avgp2quar = np.array(fourth_order.get( 'avgp2quar' ))
                self.avgp3quar = np.array(fourth_order.get( 'avgp3quar' ))
                self.avgx1sqp1sq = np.array(fourth_order.get( 'avgx1sqp1sq' ))
                self.avgx2sqp2sq = np.array(fourth_order.get( 'avgx2sqp2sq' ))
                self.avgx3sqp3sq = np.array(fourth_order.get( 'avgx3sqp3sq' ))


    def write(self, file, order = None):
        H5File.__init__(self, file)
        if order == None:
            order = self.__order

        if not self.is_h5_file():
            print('Error:\tFile name does not have an hdf5 ending!')
            sys.exit(1) 

        # TODO: write order of moments and if crossterms into attributes!!!

        h5f = h5py.File(file, "w")

        h5f.attrs['mom_order'] = order

        dset_zeta_array = h5f.create_dataset( self.zeta_array_name, data = self.zeta_array )
        dset_time_array = h5f.create_dataset( self.time_array_name, data = self.time_array )

        dset_charge = h5f.create_dataset(  "charge", data = self.charge )
        
        firor = h5f.create_group("first_order")

        dset_avgx1 = firor.create_dataset(  "avgx1", data = self.avgx1 )
        dset_avgx2 = firor.create_dataset(  "avgx2", data = self.avgx2 )
        dset_avgx3 = firor.create_dataset(  "avgx3", data = self.avgx3 )
        dset_avgp1 = firor.create_dataset(  "avgp1", data = self.avgp1 )
        dset_avgp2 = firor.create_dataset(  "avgp2", data = self.avgp2 )
        dset_avgp3 = firor.create_dataset(  "avgp3", data = self.avgp3 )

        if order>1:
            secor = h5f.create_group("second_order")
            dset_avgx1sq = secor.create_dataset(  "avgx1sq", data = self.avgx1sq )
            dset_avgx2sq = secor.create_dataset(  "avgx2sq", data = self.avgx2sq )
            dset_avgx3sq = secor.create_dataset(  "avgx3sq", data = self.avgx3sq )
            dset_avgp1sq = secor.create_dataset(  "avgp1sq", data = self.avgp1sq )
            dset_avgp2sq = secor.create_dataset(  "avgp2sq", data = self.avgp2sq )
            dset_avgp3sq = secor.create_dataset(  "avgp3sq", data = self.avgp3sq )
            dset_avgx1p1 = secor.create_dataset(  "avgx1p1", data = self.avgx1p1 )
            dset_avgx2p2 = secor.create_dataset(  "avgx2p2", data = self.avgx2p2 )
            dset_avgx3p3 = secor.create_dataset(  "avgx3p3", data = self.avgx3p3 )

            # crossterms:
            dset_avgx1x2 = secor.create_dataset(  "avgx1x2", data = self.avgx1x2 )
            dset_avgx3x1 = secor.create_dataset(  "avgx3x1", data = self.avgx3x1 )
            dset_avgx2x3 = secor.create_dataset(  "avgx2x3", data = self.avgx2x3 )
            dset_avgp1p2 = secor.create_dataset(  "avgp1p2", data = self.avgp1p2 )
            dset_avgp3p1 = secor.create_dataset(  "avgp3p1", data = self.avgp3p1 )
            dset_avgp2p3 = secor.create_dataset(  "avgp2p3", data = self.avgp2p3 )
            dset_avgx1p2 = secor.create_dataset(  "avgx1p2", data = self.avgx1p2 )
            dset_avgx1p3 = secor.create_dataset(  "avgx1p3", data = self.avgx1p3 )
            dset_avgx2p1 = secor.create_dataset(  "avgx2p1", data = self.avgx2p1 )
            dset_avgx2p3 = secor.create_dataset(  "avgx2p3", data = self.avgx2p3 )
            dset_avgx3p1 = secor.create_dataset(  "avgx3p1", data = self.avgx3p1 )
            dset_avgx3p2 = secor.create_dataset(  "avgx3p2", data = self.avgx3p2 )

        if order>2:
            thior = h5f.create_group("third_order")
            dset_avgx1cube = thior.create_dataset(  "avgx1cube", data = self.avgx1cube )
            dset_avgx2cube = thior.create_dataset(  "avgx2cube", data = self.avgx2cube )
            dset_avgx3cube = thior.create_dataset(  "avgx3cube", data = self.avgx3cube )
            dset_avgp1cube = thior.create_dataset(  "avgp1cube", data = self.avgp1cube )
            dset_avgp2cube = thior.create_dataset(  "avgp2cube", data = self.avgp2cube )
            dset_avgp3cube = thior.create_dataset(  "avgp3cube", data = self.avgp3cube )
            dset_avgx1sqp1 = thior.create_dataset(  "avgx1sqp1", data = self.avgx1sqp1 )
            dset_avgx2sqp2 = thior.create_dataset(  "avgx2sqp2", data = self.avgx2sqp2 )
            dset_avgx3sqp3 = thior.create_dataset(  "avgx3sqp3", data = self.avgx3sqp3 )
            dset_avgx1p1sq = thior.create_dataset(  "avgx1p1sq", data = self.avgx1p1sq )
            dset_avgx2p2sq = thior.create_dataset(  "avgx2p2sq", data = self.avgx2p2sq )
            dset_avgx3p3sq = thior.create_dataset(  "avgx3p3sq", data = self.avgx3p3sq )

        if order>3:
            fouor = h5f.create_group("fourth_order")
            dset_avgx1quar = fouor.create_dataset(  "avgx1quar", data = self.avgx1quar )
            dset_avgx2quar = fouor.create_dataset(  "avgx2quar", data = self.avgx2quar )
            dset_avgx3quar = fouor.create_dataset(  "avgx3quar", data = self.avgx3quar )
            dset_avgp1quar = fouor.create_dataset(  "avgp1quar", data = self.avgp1quar )
            dset_avgp2quar = fouor.create_dataset(  "avgp2quar", data = self.avgp2quar )
            dset_avgp3quar = fouor.create_dataset(  "avgp3quar", data = self.avgp3quar )
            dset_avgx1sqp1sq = fouor.create_dataset(  "avgx1sqp1sq", data = self.avgx1sqp1sq )
            dset_avgx2sqp2sq = fouor.create_dataset(  "avgx2sqp2sq", data = self.avgx2sqp2sq )
            dset_avgx3sqp3sq = fouor.create_dataset(  "avgx3sqp3sq", data = self.avgx3sqp3sq )

        h5f.close()

    def truncate_zeta_region(self, zeta_min, zeta_max, order = None, crossterms = False):

        if order == None:
            order = self.__order

        idx = np.nonzero(np.logical_and( self.zeta_array >= zeta_min, self.zeta_array <= zeta_max ))[0]

        self.zeta_array = self.zeta_array[idx]
        self.charge = self.charge[:,idx]

        if order>0:
            self.avgx1 = self.avgx1[:,idx]
            self.avgx2 = self.avgx2[:,idx]
            self.avgx3 = self.avgx3[:,idx]
            self.avgp1 = self.avgp1[:,idx]
            self.avgp2 = self.avgp2[:,idx]
            self.avgp3 = self.avgp3[:,idx]

        if order>1:
            self.avgx1sq = self.avgx1sq[:,idx]
            self.avgx2sq = self.avgx2sq[:,idx]
            self.avgx3sq = self.avgx3sq[:,idx]
            self.avgp1sq = self.avgp1sq[:,idx]
            self.avgp2sq = self.avgp2sq[:,idx]
            self.avgp3sq = self.avgp3sq[:,idx]
            self.avgx1p1 = self.avgx1p1[:,idx]
            self.avgx2p2 = self.avgx2p2[:,idx]
            self.avgx3p3 = self.avgx3p3[:,idx]

            if crossterms:
                self.avgx1x2 = self.avgx1x2[:,idx]
                self.avgx3x1 = self.avgx1x2[:,idx]
                self.avgx2x3 = self.avgx1x2[:,idx]
                self.avgp1p2 = self.avgp1p2[:,idx]
                self.avgp3p1 = self.avgp3p1[:,idx]
                self.avgp2p3 = self.avgp2p3[:,idx]
                self.avgx1p2 = self.avgx1p2[:,idx]
                self.avgx1p3 = self.avgx1p3[:,idx]
                self.avgx2p1 = self.avgx2p1[:,idx]
                self.avgx2p3 = self.avgx2p3[:,idx]
                self.avgx3p1 = self.avgx3p1[:,idx]
                self.avgx3p2 = self.avgx3p2[:,idx]

        if order>2:
            self.avgx1cube = self.avgx1cube[:,idx]
            self.avgx2cube = self.avgx2cube[:,idx]
            self.avgx3cube = self.avgx3cube[:,idx]
            self.avgp1cube = self.avgp1cube[:,idx]
            self.avgp2cube = self.avgp2cube[:,idx]
            self.avgp3cube = self.avgp3cube[:,idx]
            self.avgx1sqp1 = self.avgx1sqp1[:,idx]
            self.avgx2sqp2 = self.avgx2sqp2[:,idx]
            self.avgx3sqp3 = self.avgx3sqp3[:,idx]
            self.avgx1p1sq = self.avgx1p1sq[:,idx]
            self.avgx2p2sq = self.avgx2p2sq[:,idx]
            self.avgx3p3sq = self.avgx3p3sq[:,idx]

        if order>3:
            self.avgx1quar = self.avgx1quar[:,idx]
            self.avgx2quar = self.avgx2quar[:,idx]
            self.avgx3quar = self.avgx3quar[:,idx]
            self.avgp1quar = self.avgp1quar[:,idx]
            self.avgp2quar = self.avgp2quar[:,idx]
            self.avgp3quar = self.avgp3quar[:,idx]
            self.avgx1sqp1sq = self.avgx1sqp1sq[:,idx]
            self.avgx2sqp2sq = self.avgx2sqp2sq[:,idx]
            self.avgx3sqp3sq = self.avgx3sqp3sq[:,idx]


class H5Plot:
    def __init__(self):
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


class H5FList():
    def __init__(self, paths, h5ftype=None):

        self.paths = paths
        self.h5ftype = h5ftype
        self.flist = None

    def get(self, verbose=True):
        if not self.paths:
            print('Error: No file provided!')
            sys.exit(1)

        flist = []

        for path in self.paths:
            if os.path.isfile(path):
                file = path
                h5f = H5File(file, h5ftype=self.h5ftype)
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
                        h5f = H5File(file, h5ftype=self.h5ftype)
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
        self.flist = sorted(flist)
        return self.flist

    def get_uniques(self):
        n_time_chars = 6;
        fnames = []
        if self.flist == None:
            self.get()
        for f in self.flist:
            h5f = H5File(f)
            fnames.append(h5f.get_filename_wo_time())
        return list(set(fnames))

    def split_by_uniques(self):
        if self.flist == None:
            self.get()
        uniques = self.get_uniques()

        # initialize and append to list of lists
        lofl = [[] for i in range(len(uniques))] 
        for i in range(len(uniques)):
            for f in self.flist:
                if uniques[i] in os.path.split(f)[1]:
                    lofl[i].append(f)
        return lofl


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
