#!/usr/bin/env python3
# ps_ana.py

import os
import numpy as np
import h5py
import picdefs
import sys


def moments(array1, array2, weights, order=2, central=True, roots=False):
   
  # normalize weights
  w = weights/sum(weights);
  
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
      
  def moments(self, order=2):    
    for ibin in range(0,self.nbins):
      indices = [ i for i,x in enumerate(self.raw.x1) if 
                  ( (x >= self.edges[ibin]) & (x < self.edges[ibin+1]) ) ]
      q_sl = []
      x2_sl = []
      p2_sl = []
      for i in indices:
        q_sl.append(self.raw.q[i])
        x2_sl.append(self.raw.x2[i])
        p2_sl.append(self.raw.p2[i])
        
      print(moments(np.asarray(x2_sl), np.asarray(p2_sl), np.asarray(q_sl), order = order))
 
 
 
 
  