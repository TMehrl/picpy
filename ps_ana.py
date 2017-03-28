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
      
  def calc_moments(self, order=2, central=True):    
    
    partweight = np.zeros((self.raw.npart,), dtype=np.float32)
    
    self.charge = np.zeros((self.nbins,), dtype=np.float32)
  
    ibinpart = np.searchsorted(self.edges, self.raw.x1)     
  
    # Consider doing a searchsorted of ibinpart and accordingly sort of x and p arrays...

    for i in range(0,self.raw.npart):
      self.charge[ibinpart[i]] += self.raw.q[i]

    for i in range(0,self.raw.npart):
      partweight[i] = self.raw.q[i]/self.charge[ibinpart[i]]

    if order > 0:
      self.avgx1 = np.zeros((self.nbins, ), dtype=np.float32)
      self.avgx2 = np.zeros((self.nbins, ), dtype=np.float32)
      self.avgx3 = np.zeros((self.nbins, ), dtype=np.float32)
      self.avgp1 = np.zeros((self.nbins, ), dtype=np.float32)
      self.avgp2 = np.zeros((self.nbins, ), dtype=np.float32)
      self.avgp3 = np.zeros((self.nbins, ), dtype=np.float32)
          
      for i in range(0,self.raw.npart):
        self.avgx1[ibinpart[i]] += self.raw.x1[i] * partweight[i]   
        self.avgx2[ibinpart[i]] += self.raw.x2[i] * partweight[i]
        self.avgx3[ibinpart[i]] += self.raw.x3[i] * partweight[i]
        self.avgp1[ibinpart[i]] += self.raw.p1[i] * partweight[i]   
        self.avgp2[ibinpart[i]] += self.raw.p2[i] * partweight[i]
        self.avgp3[ibinpart[i]] += self.raw.p3[i] * partweight[i]

    if order > 1:
      self.avgx1sq = np.zeros((self.nbins, ), dtype=np.float32)
      self.avgx2sq = np.zeros((self.nbins, ), dtype=np.float32)
      self.avgx3sq = np.zeros((self.nbins, ), dtype=np.float32)
      self.avgp1sq = np.zeros((self.nbins, ), dtype=np.float32)
      self.avgp2sq = np.zeros((self.nbins, ), dtype=np.float32)
      self.avgp3sq = np.zeros((self.nbins, ), dtype=np.float32)
      self.avgx1p1 = np.zeros((self.nbins, ), dtype=np.float32)
      self.avgx2p2 = np.zeros((self.nbins, ), dtype=np.float32)
      self.avgx3p3 = np.zeros((self.nbins, ), dtype=np.float32)   
      
      for i in range(0,self.raw.npart):
        self.avgx1sq[ibinpart[i]] += self.raw.x1[i]**2 * partweight[i]   
        self.avgx2sq[ibinpart[i]] += self.raw.x2[i]**2 * partweight[i]
        self.avgx3sq[ibinpart[i]] += self.raw.x3[i]**2 * partweight[i]
        self.avgp1sq[ibinpart[i]] += self.raw.p1[i]**2 * partweight[i]   
        self.avgp2sq[ibinpart[i]] += self.raw.p2[i]**2 * partweight[i]
        self.avgp3sq[ibinpart[i]] += self.raw.p3[i]**2 * partweight[i]
        self.avgx1p1[ibinpart[i]] += self.raw.x1[i] * self.raw.p1[i] * partweight[i]
        self.avgx2p2[ibinpart[i]] += self.raw.x2[i] * self.raw.p2[i] * partweight[i]
        self.avgx3p3[ibinpart[i]] += self.raw.x3[i] * self.raw.p3[i] * partweight[i]

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
        
    if order > 2:
      print('Moments with orders > 3 not yet implemented!') 
            
    self.if_moms_calc = True
    
 
 
 
  