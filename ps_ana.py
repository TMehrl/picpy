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
  
    ibinpart_raw = np.searchsorted(self.edges, self.raw.x1)     
    self.cm_afsearchsorted_time = time.time()
    
    ibinpartsort = np.argsort(ibinpart_raw) 
    self.cm_afargsort_time = time.time()
    
    ibinpart = np.zeros((self.raw.npart, ), dtype=np.uint32)
    q = np.zeros((self.raw.npart, ), dtype=np.float32)
    x1 = np.zeros((self.raw.npart, ), dtype=np.float32)
    x2 = np.zeros((self.raw.npart, ), dtype=np.float32)
    x3 = np.zeros((self.raw.npart, ), dtype=np.float32)
    p1 = np.zeros((self.raw.npart, ), dtype=np.float32)
    p2 = np.zeros((self.raw.npart, ), dtype=np.float32)
    p3 = np.zeros((self.raw.npart, ), dtype=np.float32)
 
    self.cm_afallocsortpart_time = time.time()
    
    for i in range(0,self.raw.npart):
      ipart_sort = ibinpartsort[i]
      ibinpart[ipart_sort ] = ibinpart_raw[i]
      q[  ipart_sort ] =  self.raw.q[i]
      x1[ ipart_sort ] = self.raw.x1[i]
      x2[ ipart_sort ] = self.raw.x2[i]
      x3[ ipart_sort ] = self.raw.x3[i]
      p1[ ipart_sort ] = self.raw.p1[i]
      p2[ ipart_sort ] = self.raw.p2[i]
      p3[ ipart_sort ] = self.raw.p3[i]

    self.cm_afsortingpart_time = time.time()
    
    self.charge = np.zeros((self.nbins,), dtype=np.float32)  
    for i in range(0,self.raw.npart):
      self.charge[ibinpart[i]] += q[i]
        
    partweight = np.zeros((self.raw.npart,), dtype=np.float32)
      
    for i in range(0,self.raw.npart):
      partweight[i] = q[i]/self.charge[ibinpart[i]]
      
    self.cm_afpartweight_time = time.time()
    if order > 0:
      self.avgx1 = np.zeros((self.nbins, ), dtype=np.float32)
      self.avgx2 = np.zeros((self.nbins, ), dtype=np.float32)
      self.avgx3 = np.zeros((self.nbins, ), dtype=np.float32)
      self.avgp1 = np.zeros((self.nbins, ), dtype=np.float32)
      self.avgp2 = np.zeros((self.nbins, ), dtype=np.float32)
      self.avgp3 = np.zeros((self.nbins, ), dtype=np.float32)
      self.cm_afallocavg_time = time.time()
          
      for i in range(0,self.raw.npart):
        ibin = ibinpart[i]
        self.avgx1[ibin] += x1[i] * partweight[i]   
        self.avgx2[ibin] += x2[i] * partweight[i]
        self.avgx3[ibin] += x3[i] * partweight[i]
        self.avgp1[ibin] += p1[i] * partweight[i]   
        self.avgp2[ibin] += p2[i] * partweight[i]
        self.avgp3[ibin] += p3[i] * partweight[i]
      self.cm_afcalcavg_time = time.time()
      
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
      self.cm_afallocsqavg_time = time.time()
      
    
      part_avgx1sq = np.multiply(np.power(x1,2 ), partweight ) 
      part_avgx2sq = np.multiply(np.power(x2,2 ), partweight )
      part_avgx3sq = np.multiply(np.power(x3,2 ), partweight )
      part_avgp1sq = np.multiply(np.power(p1,2 ), partweight )
      part_avgp2sq = np.multiply(np.power(p2,2 ), partweight )
      part_avgp3sq = np.multiply(np.power(p3,2 ), partweight )
      part_avgx1p1 = np.multiply(np.multiply( x1, p1 ), partweight )
      part_avgx2p2 = np.multiply(np.multiply( x2, p2 ), partweight )
      part_avgx3p3 = np.multiply(np.multiply( x3, p3 ), partweight )
      
      for i in range(0,self.raw.npart):
        ibin = ibinpart[i]        
        self.avgx1sq[ibin] += part_avgx1sq[i]   
        self.avgx2sq[ibin] += part_avgx2sq[i]
        self.avgx3sq[ibin] += part_avgx3sq[i]
        self.avgp1sq[ibin] += part_avgp1sq[i]   
        self.avgp2sq[ibin] += part_avgp2sq[i]
        self.avgp3sq[ibin] += part_avgp3sq[i]
        self.avgx1p1[ibin] += part_avgx1p1[i]
        self.avgx2p2[ibin] += part_avgx2p2[i]
        self.avgx3p3[ibin] += part_avgx3p3[i]
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
      print('Searchsorted1:\t\t%e %s' % ((self.cm_afsearchsorted_time-self.startcm_time), 's'))
      print('Argsort:\t\t%e %s' % ((self.cm_afargsort_time-self.cm_afsearchsorted_time), 's'))
      print('Alloc sorted part arr:\t%e %s' % ((self.cm_afallocsortpart_time-self.cm_afargsort_time), 's'))
      print('Sort part arr:\t\t%e %s' % ((self.cm_afsortingpart_time-self.cm_afallocsortpart_time), 's'))
      print('Calc of partweight:\t%e %s' % ((self.cm_afpartweight_time-self.cm_afsortingpart_time), 's'))
      print('Alloc of avgs:\t\t%e %s' % ((self.cm_afallocavg_time-self.cm_afpartweight_time), 's'))
      print('Computation of avgs:\t%e %s' % ((self.cm_afcalcavg_time-self.cm_afallocavg_time), 's'))
      print('Alloc of var arrays:\t%e %s' % ((self.cm_afallocsqavg_time-self.cm_afcalcavg_time), 's'))
      print('Calc of var:\t\t%e %s' % ((self.cm_afcalcsqavg_time-self.cm_afallocsqavg_time), 's'))
      print('Calc of centr of vars:\t%e %s' % ((self.cm_afcalcsqavgcent_time-self.cm_afcalcsqavg_time), 's'))
      
    if order > 2:
      print('Moments with orders > 3 not yet implemented!') 
            
    self.if_moms_calc = True
    
 
 
 
  