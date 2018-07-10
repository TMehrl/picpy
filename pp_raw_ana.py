#!/usr/bin/env python3
# pp_raw_ana.py

import os
import numpy as np
import gc
import h5py
import pp_defs
import sys
import time


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

# Class for slice analysis
class Slices:
    def __init__(self, raw, edges=[], nbins=0, zrange=None, cellvol=1.0):

        self.edges = edges
        self.if_edges_eq_spaced = True
        self.if_moms_calc = False
        self.cellvol = cellvol

        dx0 = (raw.xmax[0] - raw.xmin[0])/raw.nx[0]
        self.raw = raw
        if nbins == 0:
            if (edges==[]) and (zrange == None):
                self.edges = np.linspace(raw.xmin[0]-dx0/2, raw.xmax[0]+dx0/2, num=(raw.nx[0]+1))
                self.if_edges_eq_spaced = True
            elif (edges==[]) and (zrange != None):
                self.edges = np.arange(start=zrange[0]-dx0/2,stop=zrange[1]+dx0/2,step=dx0,dtype=np.float32)
                self.if_edges_eq_spaced = True
        elif nbins != 0:
            if zrange==None:
                self.edges = np.linspace(raw.xmin[0]-dx0/2, raw.xmax[0]+dx0/2, num=(nbins+1))
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

    class time:
        startcm = 0.0
        cm_afcalcsqavg = 0.0
        cm_afallocsortpart = 0.0
        cm_afsortingpart = 0.0
        cm_afcalcavg = 0.0
        cm_afallocsqavg = 0.0       


    def alloc_arrays(self, nbins):

        self.avgx1 = np.zeros((nbins), dtype=np.float32)
        self.avgx2 = np.zeros((nbins), dtype=np.float32)
        self.avgx3 = np.zeros((nbins), dtype=np.float32)
        self.avgp1 = np.zeros((nbins), dtype=np.float32)
        self.avgp2 = np.zeros((nbins), dtype=np.float32)
        self.avgp3 = np.zeros((nbins), dtype=np.float32)      

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


    def calc_moments(self, order=2, central=True, crossterms=False, timings=False, reshape_method=False):

        if timings: self.time.startcm = time.time()

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

        if timings: self.time.cm_afsearchsorted = time.time()

        if timings: self.time.cm_afallocsortpart = time.time()

        # Sorting partile arrays according to the bins they are located in
        idx = np.argsort(ibinpart)
        q = q[idx]
        x1 = x1[idx]
        x2 = x2[idx]
        x3 = x3[idx]
        p1 = p1[idx]
        p2 = p2[idx]
        p3 = p3[idx]

        if timings: self.time.cm_afsortingpart = time.time()

        # Setting idx range for each bin
        i1 = np.cumsum(self.npart) - self.npart
        i2 = np.cumsum(self.npart) - 1

        if order > 0:
            for ibin in range(0,self.nbins):
                # Making sure sum of weights is not zero:
                if  self.charge[ibin] != 0.0:
                    self.avgx1[ibin] = np.ma.average(x1[i1[ibin]:i2[ibin]], weights=q[i1[ibin]:i2[ibin]])
                    self.avgx2[ibin] = np.ma.average(x2[i1[ibin]:i2[ibin]], weights=q[i1[ibin]:i2[ibin]])
                    self.avgx3[ibin] = np.ma.average(x3[i1[ibin]:i2[ibin]], weights=q[i1[ibin]:i2[ibin]])
                    self.avgp1[ibin] = np.ma.average(p1[i1[ibin]:i2[ibin]], weights=q[i1[ibin]:i2[ibin]])
                    self.avgp2[ibin] = np.ma.average(p2[i1[ibin]:i2[ibin]], weights=q[i1[ibin]:i2[ibin]])
                    self.avgp3[ibin] = np.ma.average(p3[i1[ibin]:i2[ibin]], weights=q[i1[ibin]:i2[ibin]])

            if timings: self.time.cm_afcalcavg = time.time()


        if order > 1:
            if central:
                for ibin in range(0,self.nbins):
                    x1[i1[ibin]:i2[ibin]] = x1[i1[ibin]:i2[ibin]] - self.avgx1[ibin]
                    x2[i1[ibin]:i2[ibin]] = x2[i1[ibin]:i2[ibin]] - self.avgx2[ibin]
                    x3[i1[ibin]:i2[ibin]] = x3[i1[ibin]:i2[ibin]] - self.avgx3[ibin]
                    p1[i1[ibin]:i2[ibin]] = p1[i1[ibin]:i2[ibin]] - self.avgp1[ibin]
                    p2[i1[ibin]:i2[ibin]] = p2[i1[ibin]:i2[ibin]] - self.avgp2[ibin]
                    p3[i1[ibin]:i2[ibin]] = p3[i1[ibin]:i2[ibin]] - self.avgp3[ibin]

            if timings: self.time.cm_afallocsqavg = time.time()

            for ibin in range(0,self.nbins):
                # Making sure sum of weights is not zero:
                if  self.charge[ibin] != 0.0:
                    self.avgx1sq[ibin] = np.ma.average(np.power(x1[i1[ibin]:i2[ibin]],2), weights=q[i1[ibin]:i2[ibin]])
                    self.avgx2sq[ibin] = np.ma.average(np.power(x2[i1[ibin]:i2[ibin]],2), weights=q[i1[ibin]:i2[ibin]])
                    self.avgx3sq[ibin] = np.ma.average(np.power(x3[i1[ibin]:i2[ibin]],2), weights=q[i1[ibin]:i2[ibin]])
                    self.avgp1sq[ibin] = np.ma.average(np.power(p1[i1[ibin]:i2[ibin]],2), weights=q[i1[ibin]:i2[ibin]])
                    self.avgp2sq[ibin] = np.ma.average(np.power(p2[i1[ibin]:i2[ibin]],2), weights=q[i1[ibin]:i2[ibin]])
                    self.avgp3sq[ibin] = np.ma.average(np.power(p3[i1[ibin]:i2[ibin]],2), weights=q[i1[ibin]:i2[ibin]])
       
                    self.avgx1p1[ibin] = np.ma.average(np.multiply(x1[i1[ibin]:i2[ibin]],p1[i1[ibin]:i2[ibin]]), weights=q[i1[ibin]:i2[ibin]])
                    self.avgx2p2[ibin] = np.ma.average(np.multiply(x2[i1[ibin]:i2[ibin]],p2[i1[ibin]:i2[ibin]]), weights=q[i1[ibin]:i2[ibin]])
                    self.avgx3p3[ibin] = np.ma.average(np.multiply(x3[i1[ibin]:i2[ibin]],p3[i1[ibin]:i2[ibin]]), weights=q[i1[ibin]:i2[ibin]])

            if crossterms:
                for ibin in range(0,self.nbins):
                    # Making sure sum of weights is not zero:
                    if  self.charge[ibin] != 0.0:
                        self.avgx1x2[ibin] = np.ma.average(np.multiply(x1[i1[ibin]:i2[ibin]],p2[i1[ibin]:i2[ibin]]), weights=q[i1[ibin]:i2[ibin]])
                        self.avgx3x1[ibin] = np.ma.average(np.multiply(x3[i1[ibin]:i2[ibin]],x1[i1[ibin]:i2[ibin]]), weights=q[i1[ibin]:i2[ibin]])
                        self.avgx2x3[ibin] = np.ma.average(np.multiply(x2[i1[ibin]:i2[ibin]],x3[i1[ibin]:i2[ibin]]), weights=q[i1[ibin]:i2[ibin]])

                        self.avgp1p2[ibin] = np.ma.average(np.multiply(p1[i1[ibin]:i2[ibin]],p2[i1[ibin]:i2[ibin]]), weights=q[i1[ibin]:i2[ibin]])
                        self.avgp3p1[ibin] = np.ma.average(np.multiply(p3[i1[ibin]:i2[ibin]],p1[i1[ibin]:i2[ibin]]), weights=q[i1[ibin]:i2[ibin]])
                        self.avgp2p3[ibin] = np.ma.average(np.multiply(p2[i1[ibin]:i2[ibin]],p3[i1[ibin]:i2[ibin]]), weights=q[i1[ibin]:i2[ibin]])

                        self.avgx1p2[ibin] = np.ma.average(np.multiply(x1[i1[ibin]:i2[ibin]],p2[i1[ibin]:i2[ibin]]), weights=q[i1[ibin]:i2[ibin]])
                        self.avgx1p3[ibin] = np.ma.average(np.multiply(x1[i1[ibin]:i2[ibin]],p3[i1[ibin]:i2[ibin]]), weights=q[i1[ibin]:i2[ibin]])

                        self.avgx2p1[ibin] = np.ma.average(np.multiply(x2[i1[ibin]:i2[ibin]],p1[i1[ibin]:i2[ibin]]), weights=q[i1[ibin]:i2[ibin]])
                        self.avgx2p3[ibin] = np.ma.average(np.multiply(x2[i1[ibin]:i2[ibin]],p3[i1[ibin]:i2[ibin]]), weights=q[i1[ibin]:i2[ibin]])

                        self.avgx3p1[ibin] = np.ma.average(np.multiply(x3[i1[ibin]:i2[ibin]],p1[i1[ibin]:i2[ibin]]), weights=q[i1[ibin]:i2[ibin]])
                        self.avgx3p2[ibin] = np.ma.average(np.multiply(x3[i1[ibin]:i2[ibin]],p2[i1[ibin]:i2[ibin]]), weights=q[i1[ibin]:i2[ibin]])

            if timings: self.time.cm_afcalcsqavg = time.time()

        if order > 2:
            for ibin in range(0,self.nbins):
                # Making sure sum of weights is not zero:
                if  self.charge[ibin] != 0.0:
                    self.avgx1cube[ibin] = np.ma.average(np.power( x1[i1[ibin]:i2[ibin]],3), weights=q[i1[ibin]:i2[ibin]])
                    self.avgx2cube[ibin] = np.ma.average(np.power( x2[i1[ibin]:i2[ibin]],3), weights=q[i1[ibin]:i2[ibin]])
                    self.avgx3cube[ibin] = np.ma.average(np.power( x3[i1[ibin]:i2[ibin]],3), weights=q[i1[ibin]:i2[ibin]])

                    self.avgp1cube[ibin] = np.ma.average(np.power( p1[i1[ibin]:i2[ibin]],3), weights=q[i1[ibin]:i2[ibin]])
                    self.avgp2cube[ibin] = np.ma.average(np.power( p2[i1[ibin]:i2[ibin]],3), weights=q[i1[ibin]:i2[ibin]])
                    self.avgp3cube[ibin] = np.ma.average(np.power( p3[i1[ibin]:i2[ibin]],3), weights=q[i1[ibin]:i2[ibin]])

                    self.avgx1sqp1[ibin] = np.ma.average(np.multiply(np.power(x1[i1[ibin]:i2[ibin]],2),p1[i1[ibin]:i2[ibin]]), weights=q[i1[ibin]:i2[ibin]])
                    self.avgx2sqp2[ibin] = np.ma.average(np.multiply(np.power(x2[i1[ibin]:i2[ibin]],2),p2[i1[ibin]:i2[ibin]]), weights=q[i1[ibin]:i2[ibin]])
                    self.avgx3sqp3[ibin] = np.ma.average(np.multiply(np.power(x3[i1[ibin]:i2[ibin]],2),p3[i1[ibin]:i2[ibin]]), weights=q[i1[ibin]:i2[ibin]])

                    self.avgx1p1sq[ibin] = np.ma.average(np.multiply(x1[i1[ibin]:i2[ibin]],np.power(p1[i1[ibin]:i2[ibin]],2)), weights=q[i1[ibin]:i2[ibin]])
                    self.avgx2p2sq[ibin] = np.ma.average(np.multiply(x2[i1[ibin]:i2[ibin]],np.power(p2[i1[ibin]:i2[ibin]],2)), weights=q[i1[ibin]:i2[ibin]])
                    self.avgx3p3sq[ibin] = np.ma.average(np.multiply(x3[i1[ibin]:i2[ibin]],np.power(p3[i1[ibin]:i2[ibin]],2)), weights=q[i1[ibin]:i2[ibin]])

        if timings:
            # Timing stuff
            print('--------- Timings --------- ')
            print('Total time:\t\t%0.2e %s' % ((self.time.cm_afcalcsqavg-self.time.startcm) , 's'))
            print('Searchsorted:\t\t%0.2e %s' % ((self.time.cm_afsearchsorted-self.time.startcm), 's'))
            print('Alloc sorted part arr:\t%0.2e %s' % ((self.time.cm_afallocsortpart-self.time.cm_afsearchsorted), 's'))
            print('Sort part arr:\t\t%0.2e %s' % ((self.time.cm_afsortingpart-self.time.cm_afallocsortpart), 's'))
            print('Computation of avgs:\t%0.2e %s' % ((self.time.cm_afcalcavg-self.time.cm_afsortingpart), 's'))
            print('Alloc of var arrays:\t%0.2e %s' % ((self.time.cm_afallocsqavg-self.time.cm_afcalcavg), 's'))
            print('Calc of var:\t\t%0.2e %s' % ((self.time.cm_afcalcsqavg-self.time.cm_afallocsqavg), 's'))


        self.if_moms_calc = True

# # Class to generate histograms
# class hist:
    

