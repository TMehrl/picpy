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
    charge = [None]*self.nbins
    x1p1_moms = []
    x2p2_moms = []
    x3p3_moms = []  
    for ibin in range(0, self.nbins):
      indices = [ i for i,x in enumerate(self.raw.x1) if 
                  ( (x >= self.edges[ibin]) & (x < self.edges[ibin+1]) ) ]
      npart_sl = len(indices)
      q_sl = []
      x1_sl = []
      p1_sl = []
      x2_sl = []
      p2_sl = []
      x3_sl = []
      p3_sl = []
      for i in indices:
        q_sl.append(self.raw.q[i])
        x1_sl.append(self.raw.x1[i])
        p1_sl.append(self.raw.p1[i])
        x2_sl.append(self.raw.x2[i])
        p2_sl.append(self.raw.p2[i])
        x3_sl.append(self.raw.x3[i])
        p3_sl.append(self.raw.p3[i])
      
      charge[ibin] = np.sum(q_sl)
      x1p1_moms.append( moments( np.asarray(x1_sl), 
                                      np.asarray(p1_sl), 
                                      np.asarray(q_sl), 
                                      order = order) )               
      x2p2_moms.append( moments( np.asarray(x2_sl), 
                                      np.asarray(p2_sl), 
                                      np.asarray(q_sl), 
                                      order = order) )
      x3p3_moms.append( moments( np.asarray(x3_sl), 
                                      np.asarray(p3_sl), 
                                      np.asarray(q_sl), 
                                      order = order) )                                    
    
    self.charge = np.asarray(charge)
    self.x1p1_moms = np.asarray(x1p1_moms)
    self.x2p2_moms = np.asarray(x2p2_moms)
    self.x3p3_moms = np.asarray(x3p3_moms)
    
    self.if_moms_calc = True
    
 
 
 
  