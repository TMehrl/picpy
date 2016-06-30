#!/usr/bin/env python3.5

import h5py
import numpy as np
import sys
#import getopt
import os
import csv

h5_ext_list = ['.h5','.hdf5']
ascii_ext_list = ['.csv','.txt','.ascii']
NUM_ARGS = 2

def ext_error(fname,h5_ext_list,ascii_ext_list):
	print('ERROR:\tExtension of file:', str(fname))
	print('\tis not supported!')
	print('\tAllowed hdf5-file extensions: ',list(h5_ext_list))
	print('\tAllowed csv-file extensions: ',list(ascii_ext_list))

class H5RAW3D:
  def __init__(self, h5file):
  	self.dattype = 'hdf5'
  	with h5py.File(h5file,'r') as hf:
  		self.hdf5keys = list(hf.keys())
  		self.filename = h5file
  		self.x1 = np.array(hf.get('x1'))
  		self.x2 = np.array(hf.get('x2'))
  		self.x3 = np.array(hf.get('x3'))
  		self.p1 = np.array(hf.get('p1'))
  		self.p2 = np.array(hf.get('p2'))
  		self.p3 = np.array(hf.get('p3'))
  		self.q  = np.array(hf.get('q'))
  		self.tag  = np.array(hf.get('tag'))
  		hf.close()
  def proc_tag(self):
  	return self.tag[:,0]  	
  def part_tag(self):
  	return self.tag[:,1]			  
  def printhdf5keys(self):
  	print('List of arrays in this file: \n', self.hdf5keys)
  def sort(self):
  	ind = np.lexsort((self.proc_tag(),self.part_tag()))
  	self.x1 = self.x1[ind]
  	self.x2 = self.x2[ind]
  	self.x3 = self.x3[ind]
  	self.p1 = self.p1[ind]
  	self.p2 = self.p2[ind]
  	self.p3 = self.p3[ind]
  	self.q = self.q[ind]
  	self.tag = self.tag[ind,:]	
  	
  	
class CSVRAW3D:
  def __init__(self, csvfile):
  	self.dattype = 'csv'
  	with open(csvfile, 'rt') as cf:
  		csvstrdata = csv.reader(cf, delimiter='\t', quoting=csv.QUOTE_NONNUMERIC)
  		self.data = np.array(list(csvstrdata)).astype('float')
  		self.x1 = self.data[:,0]
  		self.x2 = self.data[:,1]
  		self.x3 = self.data[:,2]
  		self.p1 = self.data[:,3]
  		self.p2 = self.data[:,4]
  		self.p3 = self.data[:,5]
  		self.q = self.data[:,6]
  		self.tag = self.data[:,7:]
  		#self.tag = [int(tag) for tag in row for row in self.data[:,7:]]
  		self.filename = csvfile
  		cf.close()
  def proc_tag(self):
  	return np.ones(len(self.tag[:,0]))
  def part_tag(self):
  	return self.tag[:,0]	
  def sort(self):
  	ind = np.lexsort((self.proc_tag(),self.part_tag()))
  	self.x1 = self.x1[ind]
  	self.x2 = self.x2[ind]
  	self.x3 = self.x3[ind]
  	self.p1 = self.p1[ind]
  	self.p2 = self.p2[ind]
  	self.p3 = self.p3[ind]
  	self.q = self.q[ind]
  	self.tag = self.tag[ind,:]	
  	
print('This script compares two files (*.h5 or *.txt)')  
    	
if len(sys.argv) != (NUM_ARGS+1):
	print('Argument List:', str(sys.argv[1:]))
	print("ERROR: This script requires two arguments!")
	print('Usage: /.raw-compare.h5 <raw-file1> <raw-file2>')
	sys.exit()

# File 1 (Iteration not needed, since always exactly 2 files)
fname1 = sys.argv[1]
fpath1, fext1 = os.path.splitext(fname1)

print('Reading', fname1)
if any(fext1 == s for s in h5_ext_list):
	raw1 = H5RAW3D(fname1)
elif any(fext1 == s for s in ascii_ext_list):
	raw1 = CSVRAW3D(fname1)
else:
	ext_error(fname1, h5_ext_list, ascii_ext_list)
	sys.exit()
print('Reading complete.') 

# File 2
fname2 = sys.argv[2]
fpath2, fext2 = os.path.splitext(fname2)

print('Reading', fname2)
if any(fext2 == s for s in h5_ext_list):
	raw2 = H5RAW3D(fname2)
elif any(fext2 == s for s in ascii_ext_list):
	raw2 = CSVRAW3D(fname2)
else:
	ext_error(fname2, h5_ext_list, ascii_ext_list)
	sys.exit()
print('Reading complete.')   


print('Before sorting:')  
print('raw1.x1=',list(raw1.x1))
print('raw1.part_tag=',list(raw1.part_tag()))
raw1.sort()
print('After sorting:')  
print('raw1.x1=',list(raw1.x1))
print('raw1.part_tag=',list(raw1.part_tag()))

# Sorting arrays first after proc_tag then after part_tag
#ind1 = np.lexsort((raw1.proc_tag(),raw1.part_tag()))

#new_x1 = [raw1.x1[i] for i in ind1]

#print('tag =',list(raw1.part_tag()))
#print('x_new =',list(new_x1))
#print('x_old =',list(raw1.x1))
