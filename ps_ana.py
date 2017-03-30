#!/usr/bin/env python3
# ps_ana.py

import os
import numpy as np
import h5py
import picdefs
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
class SLICES:
  def __init__(self, raw, edges=[]):
  
    self.stinit_time = time.time()
    
    self.raw = raw  
    if edges == []:
      dx0 = (raw.xmax[0] - raw.xmin[0])/raw.nx[0]
      self.edges = np.linspace(raw.xmin[0]-dx0/2, raw.xmax[0]+dx0/2, num=(raw.nx[0]+1))
    else:
      self.edges
    self.nbins = len(self.edges)-1  
    self.centers = self.edges[0:-1] + np.diff(self.edges)/2
    self.if_moms_calc = False
    self.fininit_time = time.time()
      
  def calc_moments(self, order=2, central=True):    
    
    self.startcm_time = time.time()    
  
    ibinpart = np.searchsorted(self.edges, self.raw.x1)
    self.npart = np.bincount(ibinpart)
    self.charge = np.bincount(ibinpart, self.raw.q)
        
    self.cm_afsearchsorted_time = time.time()
    
    self.cm_afargsort_time = time.time()
    
    if order > 0:
    
      max_npart_sl = np.max(self.npart)
    
      bincount = np.zeros((self.nbins), dtype=np.uint32)
      Q = np.zeros((self.nbins, max_npart_sl-1), dtype=np.float32)
      # Making sure that sum of weights is not zero if no particle in bin!
      Q = np.c_[ np.ones((self.nbins, 1), dtype=np.float32), Q]
      X1 = np.zeros((self.nbins, max_npart_sl), dtype=np.float32)
      X2 = np.zeros((self.nbins, max_npart_sl), dtype=np.float32)
      X3 = np.zeros((self.nbins, max_npart_sl), dtype=np.float32)
      P1 = np.zeros((self.nbins, max_npart_sl), dtype=np.float32)
      P2 = np.zeros((self.nbins, max_npart_sl), dtype=np.float32)
      P3 = np.zeros((self.nbins, max_npart_sl), dtype=np.float32)
 
      self.cm_afallocsortpart_time = time.time()
    
      for i in range(0,self.raw.npart):
        ibin = ibinpart[i]
        ipartbin = bincount[ibin]
        Q[ ibin,  ipartbin ] = self.raw.q[i]
        X1[ ibin, ipartbin ] = self.raw.x1[i]
        X2[ ibin, ipartbin ] = self.raw.x2[i]
        X3[ ibin, ipartbin ] = self.raw.x3[i]
        P1[ ibin, ipartbin ] = self.raw.p1[i]
        P2[ ibin, ipartbin ] = self.raw.p2[i]
        P3[ ibin, ipartbin ] = self.raw.p3[i]
        bincount[ibin] += 1

      self.cm_afsortingpart_time = time.time()
          
      self.avgx1 = np.average(X1, axis=1, weights=Q)
      self.avgx2 = np.average(X2, axis=1, weights=Q)
      self.avgx3 = np.average(X3, axis=1, weights=Q)
      self.avgp1 = np.average(P1, axis=1, weights=Q)
      self.avgp2 = np.average(P2, axis=1, weights=Q)
      self.avgp3 = np.average(P3, axis=1, weights=Q)

      self.cm_afcalcavg_time = time.time()
      
    if order > 1:

      self.cm_afallocsqavg_time = time.time()
      
      self.avgx1sq = np.average(np.power(X1,2), axis=1, weights=Q)
      self.avgx2sq = np.average(np.power(X2,2), axis=1, weights=Q)
      self.avgx3sq = np.average(np.power(X3,2), axis=1, weights=Q)
      self.avgp1sq = np.average(np.power(P1,2), axis=1, weights=Q)   
      self.avgp2sq = np.average(np.power(P2,2), axis=1, weights=Q)
      self.avgp3sq = np.average(np.power(P3,2), axis=1, weights=Q) 
      self.avgx1p1 = np.average(np.multiply(X1,P1), axis=1, weights=Q)
      self.avgx2p2 = np.average(np.multiply(X2,P2), axis=1, weights=Q)
      self.avgx3p3 = np.average(np.multiply(X2,P3), axis=1, weights=Q)
       
      self.cm_afcalcsqavg_time = time.time()
      if central:
        self.avgx1sq = np.subtract(self.avgx1sq, np.power(self.avgx1, 2))
        self.avgx2sq = np.subtract(self.avgx2sq, np.power(self.avgx2, 2))
        self.avgx3sq = np.subtract(self.avgx3sq, np.power(self.avgx3, 2))
        self.avgp1sq = np.subtract(self.avgp1sq, np.power(self.avgp1, 2))
        self.avgp2sq = np.subtract(self.avgp2sq, np.power(self.avgp2, 2))
        self.avgp3sq = np.subtract(self.avgp3sq, np.power(self.avgp3, 2))
        self.avgx1p1 = np.subtract(self.avgx1p1, np.multiply(self.avgx1, self.avgp1))
        self.avgx2p2 = np.subtract(self.avgx2p2, np.multiply(self.avgx2, self.avgp2))
        self.avgx3p3 = np.subtract(self.avgx3p3, np.multiply(self.avgx3, self.avgp3))
      self.cm_afcalcsqavgcent_time = time.time() 

      # Timing stuff      
      print('--------- Timings --------- ')
      print('Total time:\t\t%e %s' % ((self.cm_afcalcsqavg_time-self.stinit_time) , 's')) 
      print('Init time:\t\t%e %s' % ((self.fininit_time-self.stinit_time), 's'))
      print('Searchsorted:\t\t%e %s' % ((self.cm_afsearchsorted_time-self.startcm_time), 's'))
      print('Alloc sorted part arr:\t%e %s' % ((self.cm_afallocsortpart_time-self.cm_afsearchsorted_time), 's'))
      print('Sort part arr:\t\t%e %s' % ((self.cm_afsortingpart_time-self.cm_afallocsortpart_time), 's'))
      print('Computation of avgs:\t%e %s' % ((self.cm_afcalcavg_time-self.cm_afsortingpart_time), 's'))
      print('Alloc of var arrays:\t%e %s' % ((self.cm_afallocsqavg_time-self.cm_afcalcavg_time), 's'))
      print('Calc of var:\t\t%e %s' % ((self.cm_afcalcsqavg_time-self.cm_afallocsqavg_time), 's'))
      print('Calc of centr of vars:\t%e %s' % ((self.cm_afcalcsqavgcent_time-self.cm_afcalcsqavg_time), 's'))
      
    if order > 2:
      print('Moments with orders > 3 not yet implemented!') 
            
    self.if_moms_calc = True
    
 
 
 
  