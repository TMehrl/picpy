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
class Slices:
  def __init__(self, raw, edges=[], nbins=0):
  
    dx0 = (raw.xmax[0] - raw.xmin[0])/raw.nx[0]
    self.raw = raw  
    if (edges==[]) & (nbins==0):
      self.edges = np.linspace(raw.xmin[0]-dx0/2, raw.xmax[0]+dx0/2, num=(raw.nx[0]+1))
    elif nbins != 0:
      self.edges = np.linspace(raw.xmin[0]-dx0/2, raw.xmax[0]+dx0/2, num=(nbins+1))
    else:
      self.edges = edges
      
    self.nbins = len(self.edges)-1  
    self.centers = self.edges[0:-1] + np.diff(self.edges)/2
    self.if_moms_calc = False
      
  def calc_moments(self, order=2, central=True, crossterms=False, timings=False):    
    
    if timings: self.startcm_time = time.time()    
  
    # Assign each particle the index of the bin it is located in
    # ( subtract 1 in order to get index '0' 
    # if particle is in interval [ edges[0] edges[1] ] etc. )
    ibinpart = np.searchsorted(self.edges, self.raw.x1) - 1 
    self.npart = np.bincount(ibinpart, minlength=self.nbins)
    self.charge = np.bincount(ibinpart, weights=self.raw.q, minlength=self.nbins)
    
    if np.size(self.npart)>self.nbins | np.size(self.charge)>self.nbins:
      print('Warning: particles out of range!')
        
    if timings: self.cm_afsearchsorted_time = time.time()
    
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
      
      if timings: self.cm_afallocsortpart_time = time.time()
    
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

      if timings: self.cm_afsortingpart_time = time.time()
          
      self.avgx1 = np.average(X1, axis=1, weights=Q)
      self.avgx2 = np.average(X2, axis=1, weights=Q)
      self.avgx3 = np.average(X3, axis=1, weights=Q)
      self.avgp1 = np.average(P1, axis=1, weights=Q)
      self.avgp2 = np.average(P2, axis=1, weights=Q)
      self.avgp3 = np.average(P3, axis=1, weights=Q)

      if timings: self.cm_afcalcavg_time = time.time()
      
    if order > 1:
    
      if central:
        X1 = X1 - self.avgx1[:,None]
        X2 = X2 - self.avgx2[:,None]
        X3 = X3 - self.avgx3[:,None]
        P1 = P1 - self.avgp1[:,None]
        P2 = P2 - self.avgp2[:,None]
        P3 = P3 - self.avgp3[:,None]

      if timings: self.cm_afallocsqavg_time = time.time()
      
      self.avgx1sq = np.average(np.power(X1,2), axis=1, weights=Q)
      self.avgx2sq = np.average(np.power(X2,2), axis=1, weights=Q)
      self.avgx3sq = np.average(np.power(X3,2), axis=1, weights=Q)

      self.avgp1sq = np.average(np.power(P1,2), axis=1, weights=Q)   
      self.avgp2sq = np.average(np.power(P2,2), axis=1, weights=Q)
      self.avgp3sq = np.average(np.power(P3,2), axis=1, weights=Q) 

      self.avgx1p1 = np.average(np.multiply(X1,P1), axis=1, weights=Q)
      self.avgx2p2 = np.average(np.multiply(X2,P2), axis=1, weights=Q)
      self.avgx3p3 = np.average(np.multiply(X2,P3), axis=1, weights=Q)
       
      if crossterms:
        self.avgx1x2 = np.average(np.multiply(X1,P2), axis=1, weights=Q)
        self.avgx3x1 = np.average(np.multiply(X1,P3), axis=1, weights=Q)
        self.avgx2x3 = np.average(np.multiply(X2,P3), axis=1, weights=Q)   
 
        self.avgp1p2 = np.average(np.multiply(X1,P2), axis=1, weights=Q)
        self.avgp3p1 = np.average(np.multiply(X1,P3), axis=1, weights=Q)
        self.avgp2p3 = np.average(np.multiply(X2,P3), axis=1, weights=Q) 
                
        self.avgx1p2 = np.average(np.multiply(X1,P2), axis=1, weights=Q)
        self.avgx1p3 = np.average(np.multiply(X1,P3), axis=1, weights=Q)

        self.avgx2p1 = np.average(np.multiply(X2,P1), axis=1, weights=Q) 
        self.avgx2p3 = np.average(np.multiply(X2,P3), axis=1, weights=Q)       

        self.avgx3p1 = np.average(np.multiply(X3,P1), axis=1, weights=Q)        
        self.avgx3p2 = np.average(np.multiply(X3,P2), axis=1, weights=Q)
        
      if timings: self.cm_afcalcsqavg_time = time.time()
    
    if order > 2:
    
      self.avgx1cube = np.average(np.power( X1 ,3), axis=1, weights=Q)
      self.avgx2cube = np.average(np.power( X2 ,3), axis=1, weights=Q)
      self.avgx3cube = np.average(np.power( X3 ,3), axis=1, weights=Q)

      self.avgp1cube = np.average(np.power( P1 ,3), axis=1, weights=Q)   
      self.avgp2cube = np.average(np.power( P2 ,3), axis=1, weights=Q)
      self.avgp3cube = np.average(np.power( P3 ,3), axis=1, weights=Q) 

      self.avgx1sqp1 = np.average(np.multiply(np.power(X1,2),P1), axis=1, weights=Q)
      self.avgx2sqp2 = np.average(np.multiply(np.power(X2,2),P2), axis=1, weights=Q)
      self.avgx3sqp3 = np.average(np.multiply(np.power(X3,2),P3), axis=1, weights=Q)

      self.avgx1p1sq = np.average(np.multiply(X1,np.power(P1,2)), axis=1, weights=Q)
      self.avgx2p2sq = np.average(np.multiply(X2,np.power(P2,2)), axis=1, weights=Q)
      self.avgx3p3sq = np.average(np.multiply(X3,np.power(P3,2)), axis=1, weights=Q)
      
    
    if timings:
      # Timing stuff      
      print('--------- Timings --------- ')
      print('Total time:\t\t%e %s' % ((self.cm_afcalcsqavg_time-self.startcm_time) , 's')) 
      print('Searchsorted:\t\t%e %s' % ((self.cm_afsearchsorted_time-self.startcm_time), 's'))
      print('Alloc sorted part arr:\t%e %s' % ((self.cm_afallocsortpart_time-self.cm_afsearchsorted_time), 's'))
      print('Sort part arr:\t\t%e %s' % ((self.cm_afsortingpart_time-self.cm_afallocsortpart_time), 's'))
      print('Computation of avgs:\t%e %s' % ((self.cm_afcalcavg_time-self.cm_afsortingpart_time), 's'))
      print('Alloc of var arrays:\t%e %s' % ((self.cm_afallocsqavg_time-self.cm_afcalcavg_time), 's'))
      print('Calc of var:\t\t%e %s' % ((self.cm_afcalcsqavg_time-self.cm_afallocsqavg_time), 's'))
      
            
    self.if_moms_calc = True
    
 
 
 
  