#!/usr/bin/env python3
# ps_ana.py

import os
import numpy as np
import h5py
import picdefs
import sys


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
    self.raw = raw  
    if edges == []:
      dx0 = (raw.xmax[0] - raw.xmin[0])/raw.nx[0]
      self.edges = np.linspace(raw.xmin[0]-dx0/2, raw.xmax[0]+dx0/2, num=(raw.nx[0]+1))
    else:
      self.edges
    self.nbins = len(self.edges)-1  
    self.centers = self.edges[0:-1] + np.diff(self.edges)/2
    self.if_moms_calc = False
      
  def calc_moments(self, order=2):    
    self.charge = np.zeros((self.nbins,), dtype=np.float32)
    self.x1p1_moms = np.zeros((self.nbins, order+1, ), dtype=np.float32)
    self.x2p2_moms = np.zeros((self.nbins, order+1, ), dtype=np.float32)
    self.x3p3_moms = np.zeros((self.nbins, order+1, ), dtype=np.float32)  
    for ibin in range(0, self.nbins):
      indices = [ i for i, x1 in enumerate(self.raw.x1) if 
                  ( (x1 >= self.edges[ibin]) & (x1 < self.edges[ibin+1]) ) ]
      npart_sl = len(indices)
      q_sl = np.zeros((npart_sl, ), dtype=np.float32)
      x1_sl = np.zeros((npart_sl, ), dtype=np.float32)
      p1_sl = np.zeros((npart_sl, ), dtype=np.float32)
      x2_sl = np.zeros((npart_sl, ), dtype=np.float32)
      p2_sl = np.zeros((npart_sl, ), dtype=np.float32)
      x3_sl = np.zeros((npart_sl, ), dtype=np.float32)
      p3_sl = np.zeros((npart_sl, ), dtype=np.float32)
      
      ipart_sl = 0
      for i in indices:
        q_sl[ipart_sl] = self.raw.q[i]
        x1_sl[ipart_sl] = self.raw.x1[i]
        p1_sl[ipart_sl] = self.raw.p1[i]
        x2_sl[ipart_sl] = self.raw.x2[i]
        p2_sl[ipart_sl] = self.raw.p2[i]
        x3_sl[ipart_sl] = self.raw.x3[i]
        p3_sl[ipart_sl] = self.raw.p3[i]
        ipart_sl += 1
      
      self.charge[ibin] = np.sum(q_sl)
      self.x1p1_moms[ibin,:] = moments( x1_sl, 
                                        p1_sl, 
                                        q_sl, 
                                        order = order)               
      self.x2p2_moms[ibin,:] = moments( x2_sl, 
                                        p2_sl, 
                                        q_sl, 
                                        order = order)
      self.x3p3_moms[ibin,:] = moments( x3_sl, 
                                        p3_sl, 
                                        q_sl, 
                                        order = order)                                     
    
    self.if_moms_calc = True
    
 
 
 
  