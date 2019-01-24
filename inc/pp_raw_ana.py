#!/usr/bin/env python3
# pp_raw_ana.py

import os
import sys
import numpy as np
from numba import jit
import gc
import h5py
import time
import pp_defs
from pp_h5dat import HiRAW

def moments(array1, array2, weights, order=2, central=True, roots=False):

    # normalize weights
    w = weights/np.sum(weights);

    # compute means
    mean = [  np.dot( array1, w),
              np.dot( array2, w) ]

    if central:
        a1 = array1 - mean[0]
        a2 = array2 - mean[1]
    else:
        a1 = array1
        a2 = array2

    if order == 1:
        mom = mean
    elif order == 2:
        mom = [ np.dot( np.multiply(a1, a1), w ),
                np.dot( np.multiply(a1, a2), w ),
                np.dot( np.multiply(a2, a2), w ) ]
    elif order == 3:
        mom = [ np.dot( np.multiply( np.multiply(a1, a1), a1) , w ),
                np.dot( np.multiply( np.multiply(a1, a1), a2) , w ),
                np.dot( np.multiply( np.multiply(a1, a2), a2) , w ),
                np.dot( np.multiply( np.multiply(a2, a2), a2) , w ) ]
    else:
        print('Moments with orders > 3 not yet implemented!')

    if roots:
        return np.power(mom,1/order)
    else:
        return mom


@jit(nopython=True) # Set "nopython" mode for best performance
def sl_mom1(psv1,weights,i1,i2):
    Nbins = i1.size
    mom1 = np.zeros((Nbins), dtype=np.float32)
    for ibin in range(0,Nbins):
        i1b = i1[ibin]
        i2b = i2[ibin]
        weight_sum = np.sum(weights[i1b:i2b])
        if weight_sum > 0:
            mom1[ibin] = np.dot(psv1[i1b:i2b], weights[i1b:i2b])/weight_sum
    return mom1

@jit(nopython=True) # Set "nopython" mode for best performance
def sl_mom2(psv1,psv2,weights,i1,i2):
    Nbins = i1.size
    mom2 = np.zeros((Nbins), dtype=np.float32)
    for ibin in range(0,Nbins):
        i1b = i1[ibin]
        i2b = i2[ibin]
        weight_sum = np.sum(weights[i1b:i2b])
        if weight_sum > 0:
            mom2[ibin] = np.dot( np.multiply(psv1[i1b:i2b],psv2[i1b:i2b]), weights[i1b:i2b])/weight_sum
    return mom2


@jit(nopython=True) # Set "nopython" mode for best performance
def sl_mom3(psv1,psv2,psv3,weights,i1,i2):
    Nbins = i1.size
    mom3 = np.zeros((Nbins), dtype=np.float32)
    for ibin in range(0,Nbins):
        i1b = i1[ibin]
        i2b = i2[ibin]
        weight_sum = np.sum(weights[i1b:i2b])
        if weight_sum > 0:
            mom3[ibin] = np.dot( np.multiply(np.multiply(psv1[i1b:i2b],psv2[i1b:i2b]),psv3[i1b:i2b]), weights[i1b:i2b])/weight_sum
    return mom3


@jit(nopython=True) # Set "nopython" mode for best performance
def sl_mom4(psv1,psv2,psv3,psv4,weights,i1,i2):
    Nbins = i1.size
    mom4 = np.zeros((Nbins), dtype=np.float32)
    for ibin in range(0,Nbins):
        i1b = i1[ibin]
        i2b = i2[ibin]
        weight_sum = np.sum(weights[i1b:i2b])
        if weight_sum > 0:
            mom4[ibin] = np.dot( np.multiply(np.multiply(psv1[i1b:i2b],psv2[i1b:i2b]),np.multiply(psv3[i1b:i2b],psv4[i1b:i2b])), weights[i1b:i2b])/weight_sum
    return mom4


@jit(nopython=True) # Set "nopython" mode for best performance
def sort_part_bin(x1,x2,x3,p1,p2,p3,q,ibinpart):
        # Sorting partile arrays according to the bins they are located in
        idx = np.argsort(ibinpart)
        return (x1[idx], x2[idx], x3[idx], p1[idx], p2[idx], p3[idx], q[idx])



# Class for slice analysis
class Slices:
    def __init__(self, raw, edges=[], nbins=0, zrange=None, cellvol=1.0):

        self.max_order = 4
        self.edges = edges
        self.if_edges_eq_spaced = True
        self.if_moms_calc = False
        self.cellvol = cellvol

        self.raw = raw
        dx0 = self.raw.get_dx(0)        
        if nbins == 0:
            if (edges==[]) and (zrange == None):
                self.edges = np.linspace(raw.get_xmin(0)-dx0/2, raw.get_xmax(0)+dx0/2, num=(raw.get_nx(0)+1))
                self.if_edges_eq_spaced = True
            elif (edges==[]) and (zrange != None):
                self.edges = np.arange(start=zrange[0]-dx0/2,stop=zrange[1]+dx0/2,step=dx0,dtype=np.float32)
                self.if_edges_eq_spaced = True
        elif nbins != 0:
            if zrange==None:
                self.edges = np.linspace(raw.get_xmin(0)-dx0/2, raw.get_xmax(0)+dx0/2, num=(nbins+1))
                self.if_edges_eq_spaced = True
            else:
                dx0_step = (zrange[1] - zrange[0])/nbins
                self.edges = np.linspace(zrange[0]-dx0_step/2, zrange[1]+dx0_step/2, num=(nbins+1))
                self.if_edges_eq_spaced = True
        else:
            self.edges = edges
            self.if_edges_eq_spaced = False

        if not np.all(np.diff(self.edges) > 0.0):
            print('Error: Bin edges are not monotonically increasing!')
            sys.exit()   

        self.dx0 = self.edges[1] - self.edges[0]
        self.nbins = len(self.edges)-1
        self.centers = self.edges[0:-1] + np.diff(self.edges)/2 

        self.alloc_arrays(nbins=self.nbins)

    class timings:
        def __init__(self, max_order):        
            self.startcm = 0.0
            self.cm_afsearchsorted = 0.0            
            self.cm_afsortingpart = 0.0    
            self.endcm = 0.0
            self.avg0 = 0.0
            self.avg1 = 0.0
            self.avg2 = 0.0
            self.avg3 = 0.0
            self.avg4 = 0.0

    def alloc_arrays(self, nbins):

        # First order
        self.avgx1 = np.zeros((nbins), dtype=np.float32)
        self.avgx2 = np.zeros((nbins), dtype=np.float32)
        self.avgx3 = np.zeros((nbins), dtype=np.float32)
        self.avgp1 = np.zeros((nbins), dtype=np.float32)
        self.avgp2 = np.zeros((nbins), dtype=np.float32)
        self.avgp3 = np.zeros((nbins), dtype=np.float32)      

        # Second order
        self.avgx1sq = np.zeros((nbins), dtype=np.float32)
        self.avgx2sq = np.zeros((nbins), dtype=np.float32)
        self.avgx3sq = np.zeros((nbins), dtype=np.float32)
        self.avgp1sq = np.zeros((nbins), dtype=np.float32)
        self.avgp2sq = np.zeros((nbins), dtype=np.float32)
        self.avgp3sq = np.zeros((nbins), dtype=np.float32) 

        self.avgx1p1 = np.zeros((nbins), dtype=np.float32)
        self.avgx2p2 = np.zeros((nbins), dtype=np.float32)
        self.avgx3p3 = np.zeros((nbins), dtype=np.float32)

        self.avgx1x2 = np.zeros((nbins), dtype=np.float32)
        self.avgx3x1 = np.zeros((nbins), dtype=np.float32)
        self.avgx2x3 = np.zeros((nbins), dtype=np.float32)

        self.avgp1p2 = np.zeros((nbins), dtype=np.float32)
        self.avgp3p1 = np.zeros((nbins), dtype=np.float32)
        self.avgp2p3 = np.zeros((nbins), dtype=np.float32)

        self.avgx1p2 = np.zeros((nbins), dtype=np.float32)
        self.avgx1p3 = np.zeros((nbins), dtype=np.float32)

        self.avgx2p1 = np.zeros((nbins), dtype=np.float32)
        self.avgx2p3 = np.zeros((nbins), dtype=np.float32)

        self.avgx3p1 = np.zeros((nbins), dtype=np.float32)
        self.avgx3p2 = np.zeros((nbins), dtype=np.float32)      

        # Third order
        self.avgx1cube = np.zeros((nbins), dtype=np.float32)
        self.avgx2cube = np.zeros((nbins), dtype=np.float32)
        self.avgx3cube = np.zeros((nbins), dtype=np.float32)

        self.avgp1cube = np.zeros((nbins), dtype=np.float32)
        self.avgp2cube = np.zeros((nbins), dtype=np.float32)
        self.avgp3cube = np.zeros((nbins), dtype=np.float32)

        self.avgx1sqp1 = np.zeros((nbins), dtype=np.float32)
        self.avgx2sqp2 = np.zeros((nbins), dtype=np.float32)
        self.avgx3sqp3 = np.zeros((nbins), dtype=np.float32)

        self.avgx1p1sq = np.zeros((nbins), dtype=np.float32)
        self.avgx2p2sq = np.zeros((nbins), dtype=np.float32)
        self.avgx3p3sq = np.zeros((nbins), dtype=np.float32) 

        # Fourth order
        self.avgx1quar = np.zeros((nbins), dtype=np.float32)
        self.avgx2quar = np.zeros((nbins), dtype=np.float32)
        self.avgx3quar = np.zeros((nbins), dtype=np.float32)

        self.avgp1quar = np.zeros((nbins), dtype=np.float32)
        self.avgp2quar = np.zeros((nbins), dtype=np.float32)
        self.avgp3quar = np.zeros((nbins), dtype=np.float32)

        self.avgx1sqp1sq = np.zeros((nbins), dtype=np.float32)
        self.avgx2sqp2sq = np.zeros((nbins), dtype=np.float32)
        self.avgx3sqp3sq = np.zeros((nbins), dtype=np.float32)


    def calc_moments(self, order=2, central=True, crossterms=False, showtimings=False, reshape_method=False):

        if showtimings: 
            timings = self.timings(self.max_order)
            timings.startcm = time.time()

        # Select subset of particles which are in range
        idx_part_in_range = np.logical_and(self.raw.x1 > self.edges[0], 
                                           self.raw.x1 <= self.edges[-1])
        x1 = self.raw.x1[idx_part_in_range]
        x2 = self.raw.x2[idx_part_in_range]
        x3 = self.raw.x3[idx_part_in_range]
        p1 = self.raw.p1[idx_part_in_range]
        p2 = self.raw.p2[idx_part_in_range]
        p3 = self.raw.p3[idx_part_in_range]
        q = self.raw.q[idx_part_in_range] * self.cellvol

        # Assign each particle the index of the bin it is located in
        if self.if_edges_eq_spaced:
            # If edges are equally spaced:
            ibinpart = np.ndarray.astype((x1 - self.edges[0])/self.dx0, dtype=np.uint32)
        else:    
            # If edges are not equally spaced: use np.searchsorted
            # ( subtract 1 in order to get index '0'
            # if particle is in interval [ edges[0] edges[1] ] etc. )
            ibinpart = np.searchsorted(self.edges, x1) - 1

        self.npart = np.bincount(ibinpart, minlength=self.nbins)
        self.charge = np.bincount(ibinpart, weights=q, minlength=self.nbins)

        if np.size(self.npart)>self.nbins | np.size(self.charge)>self.nbins:
            print('Warning: particles out of range!')

        if showtimings: 
            timings.cm_afsearchsorted = time.time()

        # Sorting partile arrays according to the bins they are located in
        idx = np.argsort(ibinpart)
        q = q[idx]
        x1 = x1[idx]
        x2 = x2[idx]
        x3 = x3[idx]
        p1 = p1[idx]
        p2 = p2[idx]
        p3 = p3[idx]

        # This seems to be slower:
        # Bottleneck is moving of data, not computing!
        # x1, x2, x3, p1, p2, p3, q = sort_part_bin(x1,x2,x3,p1,p2,p3,q,ibinpart)

        if showtimings: 
            timings.cm_afsortingpart = time.time()

        # Setting idx range for each bin
        i1 = np.cumsum(self.npart) - self.npart
        i2 = np.cumsum(self.npart) - 1

        if showtimings: 
            timings.avg0 = time.time()

        if order > 0:
            self.avgx1 = sl_mom1(x1,q,i1,i2)
            self.avgx2 = sl_mom1(x2,q,i1,i2)
            self.avgx3 = sl_mom1(x3,q,i1,i2)
            self.avgp1 = sl_mom1(p1,q,i1,i2)
            self.avgp2 = sl_mom1(p2,q,i1,i2)
            self.avgp3 = sl_mom1(p3,q,i1,i2)

        if showtimings: 
            timings.avg1 = time.time()

        if order > 1:
            if central:
                for ibin in range(0,self.nbins):
                    x1[i1[ibin]:i2[ibin]] = x1[i1[ibin]:i2[ibin]] - self.avgx1[ibin]
                    x2[i1[ibin]:i2[ibin]] = x2[i1[ibin]:i2[ibin]] - self.avgx2[ibin]
                    x3[i1[ibin]:i2[ibin]] = x3[i1[ibin]:i2[ibin]] - self.avgx3[ibin]
                    p1[i1[ibin]:i2[ibin]] = p1[i1[ibin]:i2[ibin]] - self.avgp1[ibin]
                    p2[i1[ibin]:i2[ibin]] = p2[i1[ibin]:i2[ibin]] - self.avgp2[ibin]
                    p3[i1[ibin]:i2[ibin]] = p3[i1[ibin]:i2[ibin]] - self.avgp3[ibin]

            self.avgx1sq = sl_mom2(x1,x1,q,i1,i2)
            self.avgx2sq = sl_mom2(x2,x2,q,i1,i2)
            self.avgx2sq = sl_mom2(x2,x2,q,i1,i2)
            self.avgx3sq = sl_mom2(x3,x3,q,i1,i2)
            self.avgp1sq = sl_mom2(p1,p1,q,i1,i2)
            self.avgp2sq = sl_mom2(p2,p2,q,i1,i2)
            self.avgp3sq = sl_mom2(p3,p3,q,i1,i2)

            self.avgx1p1 = sl_mom2(x1,p1,q,i1,i2)
            self.avgx2p2 = sl_mom2(x2,p2,q,i1,i2)
            self.avgx3p3 = sl_mom2(x3,p3,q,i1,i2)

            if crossterms:
                self.avgx1x2 = sl_mom2(x1,x2,q,i1,i2)
                self.avgx3x1 = sl_mom2(x3,x1,q,i1,i2)
                self.avgx2x3 = sl_mom2(x2,x3,q,i1,i2)

                self.avgp1p2 = sl_mom2(p1,p2,q,i1,i2)
                self.avgp3p1 = sl_mom2(p3,p1,q,i1,i2)
                self.avgp2p3 = sl_mom2(p2,p3,q,i1,i2)

                self.avgx1p2 = sl_mom2(x1,p2,q,i1,i2)
                self.avgx1p3 = sl_mom2(x1,p3,q,i1,i2)

                self.avgx2p1 = sl_mom2(x2,p1,q,i1,i2)
                self.avgx2p3 = sl_mom2(x2,p3,q,i1,i2)

                self.avgx3p1 = sl_mom2(x3,p1,q,i1,i2)
                self.avgx3p2 = sl_mom2(x3,p2,q,i1,i2)


        if showtimings: 
            timings.avg2 = time.time()

        if order > 2:
            self.avgx1cube = sl_mom3(x1,x1,x1,q,i1,i2)
            self.avgx2cube = sl_mom3(x2,x2,x2,q,i1,i2)
            self.avgx3cube = sl_mom3(x3,x3,x3,q,i1,i2)

            self.avgp1cube = sl_mom3(p1,p1,p1,q,i1,i2)
            self.avgp2cube = sl_mom3(p2,p2,p2,q,i1,i2)
            self.avgp3cube = sl_mom3(p3,p3,p3,q,i1,i2)

            self.avgx1sqp1 = sl_mom3(x1,x1,p1,q,i1,i2)
            self.avgx2sqp2 = sl_mom3(x2,x2,p2,q,i1,i2)
            self.avgx3sqp3 = sl_mom3(x3,x3,p3,q,i1,i2)

            self.avgx1p1sq = sl_mom3(x1,p1,p1,q,i1,i2)
            self.avgx2p2sq = sl_mom3(x2,p2,p2,q,i1,i2)
            self.avgx3p3sq = sl_mom3(x3,p3,p3,q,i1,i2)           


        if showtimings: 
            timings.avg3 = time.time()         

        if order > 3:
            self.avgx1quar = sl_mom4(x1,x1,x1,x1,q,i1,i2) 
            self.avgx2quar = sl_mom4(x2,x2,x2,x2,q,i1,i2)
            self.avgx3quar = sl_mom4(x3,x3,x3,x3,q,i1,i2)

            self.avgp1quar = sl_mom4(p1,p1,p1,p1,q,i1,i2)
            self.avgp2quar = sl_mom4(p2,p2,p2,p2,q,i1,i2)
            self.avgp3quar = sl_mom4(p3,p3,p3,p3,q,i1,i2)

            self.avgx1sqp1sq = sl_mom4(x1,x1,p1,p1,q,i1,i2)
            self.avgx2sqp2sq = sl_mom4(x2,x2,p2,p2,q,i1,i2)
            self.avgx3sqp3sq = sl_mom4(x3,x3,p3,p3,q,i1,i2)
            
        if showtimings: 
            timings.avg4 = time.time()
     

        if showtimings:
            timings.endcm = time.time() 
            # Timing stuff
            print('--------- Timings --------- ')
            print('Total time:\t\t%0.2e %s' % ((timings.endcm-timings.startcm) , 's'))
            print('Searchsorted:\t\t%0.2e %s' % ((timings.cm_afsearchsorted-timings.startcm), 's'))
            print('Sort part arr:\t\t%0.2e %s' % ((timings.cm_afsortingpart-timings.cm_afsearchsorted), 's'))
            if order > 0:
                print('Calc. 1st order moms:\t%0.2e %s' % ((timings.avg1-timings.avg0), 's'))
            if order > 1:
                print('Calc. 2nd order moms:\t%0.2e %s' % ((timings.avg2-timings.avg1), 's'))
            if order > 2:
                print('Calc. 3rd order moms:\t%0.2e %s' % ((timings.avg3-timings.avg2), 's'))
            if order > 3:
                print('Calc. 4th order moms:\t%0.2e %s' % ((timings.avg4-timings.avg3), 's'))
        self.if_moms_calc = True




class TagSelect:
    def __init__(self):
        self.raw = None
        self.ipart_sel = np.empty(0,dtype=np.int32)
        self.iproc_sel = np.empty(0,dtype=np.int32)

    def __conditional_overwrite_raw(self,raw,overwrite=False):
        print('Raw already imported!')
        if overwrite:
            print('Overwriting raw object...')
            self.raw = raw
        else:
            print('Set "overwrite=True" to overwrite raw object!')

    def import_raw_from_h5_file(self,file,verbose=True,overwrite=False):
        raw = HiRAW(file)
        raw.read_data(verbose=verbose)
        if self.raw != None:
            self.__conditional_overwrit_raw(raw,overwrite=overwrite)
        else:
            self.raw = raw           

    def import_raw_struct(self,raw,overwrite=False):
        if self.raw != None:
            self.__conditional_overwrit_raw(raw, overwrite=overwrite)
        else:
            self.raw = raw            

    # def __no_parts_selected():

    def select_part(self, selection_criterion='', overwrite=False, combine=False):
        if self.raw == None:
            print('Error: raw object not yet imported!')
        else:
            if selection_criterion == '':
                ipart_sel = self.raw.ipart
                iproc_sel = self.raw.iproc
            else:                
                x1 = self.raw.x1
                x2 = self.raw.x2
                x3 = self.raw.x3
                p1 = self.raw.p1
                p2 = self.raw.p2
                p3 = self.raw.p3

                z = self.raw.x1
                x = self.raw.x2
                y = self.raw.x3
                pz = self.raw.p1
                px = self.raw.p2
                py = self.raw.p3

                q = self.raw.q

                idx_sel = np.nonzero( eval(selection_criterion) )
                ipart_sel = self.raw.ipart[idx_sel]
                iproc_sel = self.raw.iproc[idx_sel]

            if len(self.iproc_sel) == 0 and len(self.ipart_sel) == 0:
                self.iproc_sel = iproc_sel
                self.ipart_sel = ipart_sel
            else:       
                self.__cond_owrt_cmb_tags(  iproc_sel=iproc_sel,\
                                            ipart_sel=ipart_sel,\
                                            overwrite=overwrite,\
                                            combine=combine)

    def write_taglist(self, tlname):
        if len(self.iproc_sel) == 0 or len(self.ipart_sel) == 0:
            print('Warning: No tags selected!') 
        else:
            # combine to 2d matrix
            M = np.vstack((self.iproc_sel, self.ipart_sel)).T
            np.savetxt(tlname, M, fmt='%d', delimiter=',')

    def _uniquify_tags(self):

        # Removing redundant occurences
        ipartmax = np.amax(self.ipart_sel)

        iunique = self.iproc_sel * (ipartmax+1) + self.ipart_sel

        _, idx_uniques = np.unique(iunique,return_index=True)

        self.iproc_sel = self.iproc_sel[idx_uniques]
        self.ipart_sel = self.ipart_sel[idx_uniques]

    def __cond_owrt_cmb_tags(self,iproc_sel,ipart_sel,overwrite=False,combine=False):
        print('Tags already selected!')
        
        if not overwrite and not combine:
            print('Set "overwrite=True" to overwrite tags, set "combine=True" to combine tags.')

        if overwrite:
            print('Overwriting tags...')
            self.iproc_sel = iproc_sel
            self.ipart_sel = ipart_sel

        if combine:
            print('Combining tags...')
            np.append(self.iproc_sel, iproc_sel)
            np.append(self.ipart_sel, ipart_sel)
            self._uniquify_tags()

    def read_taglist(self, tlname,overwrite=False,combine=False):
        M = np.genfromtxt(tlname, delimiter=',')
        iproc_sel = M[:,0]
        ipart_sel = M[:,1]
        
        if len(self.iproc_sel) == 0 and len(self.ipart_sel) == 0:
            self.iproc_sel = iproc_sel
            self.ipart_sel = ipart_sel
        else:       
            self.__cond_owrt_cmb_tags(  iproc_sel=iproc_sel,\
                                        ipart_sel=ipart_sel,\
                                        overwrite=overwrite,\
                                        combine=combine)

    def get_ntags(self):
        return len(self.ipart_sel)

    def get_matching_idx(self, iproc_ext_in, ipart_ext_in):
        if len(self.iproc_sel) == 0 or len(self.ipart_sel) == 0:
            print('Warning: No tags selected!') 
        else:
            ipartmax = max(np.amax(ipart_ext_in),np.amax(self.ipart_sel))

            iunique_sel = self.iproc_sel * (ipartmax+1) + self.ipart_sel
            iunique_ext = iproc_ext_in * (ipartmax+1) + ipart_ext_in

            _, idx, _ = np.intersect1d(iunique_ext,iunique_sel,assume_unique=True,return_indices=True)

            return idx


