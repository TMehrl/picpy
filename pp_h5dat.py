#!/usr/bin/env python3
# pp_h5dat.py

import os
import sys
import numpy as np
import h5py
import pp_defs


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

    def get_nx(self,dim):
        return self.nx[dim]        

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
            self.charge = np.array(hf.get( 'charge' ))
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

    def inherit_matplotlib_line_plots(self, ax):
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
