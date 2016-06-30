#!/usr/bin/env python3.5

import h5py
import numpy as np
import sys, getopt
import os
import csv

h5_ext_list = ['.h5','.hdf5']
ascii_ext_list = ['.csv','.txt','.ascii']
NUM_ARGS = 2

class H5RAW3D:
  def __init__(self, h5file):
  	self.dattype = 'hdf5'
  	with h5py.File(h5file,'r') as hf:
  		self.hdf5keys = list(hf.keys())
  		self.filename = h5file
  		self.data = np.c_[ np.array(hf.get('x1')), \
  					np.array(hf.get('x2')), \
  					np.array(hf.get('x3')), \
  					np.array(hf.get('p1')), \
  					np.array(hf.get('p2')), \
  					np.array(hf.get('p3')), \
  					np.array(hf.get('q')), \
  					np.array(hf.get('tag')) ] 
  		hf.close()
  def x1(self):
  	return self.data[:,0]
  def x2(self):
  	return self.data[:,1]	
  def x3(self):
  	return self.data[:,2]
  def p1(self):
  	return self.data[:,3]  	
  def p2(self):
  	return self.data[:,4]  
  def p3(self):
  	return self.data[:,5]  
  def q(self):
  	return self.data[:,6]  
  def tag(self):
  	return self.data[:,7:]
  def proc_tag(self):
  	return self.data[:,7]  	
  def part_tag(self):
  	return self.data[:,8]			  
  def printhdf5keys(self):
  	print('List of arrays in this file: \n', self.hdf5keys)

class CSVRAW3D:
  def __init__(self, csvfile):
  	self.dattype = 'csv'
  	with open(csvfile, 'rt') as cf:
  		csvstrdata = csv.reader(cf, delimiter='\t', quoting=csv.QUOTE_NONNUMERIC)
  		self.data = np.array(list(csvstrdata)).astype('float')
  		self.filename = csvfile
  		cf.close()
  def x1(self):
  	return self.data[:,0]
  def x2(self):
  	return self.data[:,1]	
  def x3(self):
  	return self.data[:,2]
  def p1(self):
  	return self.data[:,3]  	
  def p2(self):
  	return self.data[:,4]  
  def p3(self):
  	return self.data[:,5]  
  def q(self):
  	return self.data[:,6]  
  def tag(self):
  	return self.data[:,7:]
  def proc_tag(self):
  	return np.ones(len(self.data[:,7]))
  def part_tag(self):
  	return self.data[:,7]	
  	
  	
print('This script compares two files (*.h5 or *.txt)')  
    	
if len(sys.argv) != (NUM_ARGS+1):
	print('Argument List:', str(sys.argv[1:]))
	print("ERROR: This script requires two arguments!")
	print('Usage: /.raw-compare.h5 <raw-file1> <raw-file2>')
	sys.exit()

for arg in sys.argv[1:]:
	fname, fext = os.path.splitext(arg)
	if not( any(fext == s for s in h5_ext_list) or \
	any(fext == s for s in ascii_ext_list) ) :
		print('File:', str(arg))
		print('ERROR: File extension is not supported!')
		print('Allowed hdf5-file extensions: ',list(h5_ext_list))
		print('Allowed csv-file extensions: ',list(ascii_ext_list))
		sys.exit()

# File 1 (Iteration not needed, since always exactly 2)
fname1 = sys.argv[1]
fpath1, fext1 = os.path.splitext(fname1)

print('Reading', fname1)
if any(fext1 == s for s in h5_ext_list):
	raw1 = H5RAW3D(fname1)
elif any(fext1 == s for s in ascii_ext_list):
	raw1 = CSVRAW3D(fname1)
else:
	print('ERROR:\t\tExtension of file:', str(fname1))
	print('\t\t\tis not supported!')
	print('Allowed hdf5-file extensions: ',list(h5_ext_list))
	print('Allowed csv-file extensions: ',list(ascii_ext_list))
	sys.exit()
	
print('x1=',list(raw1.x1()))
print('tag=',list(raw1.tag()))

print('x1?: ',list(raw1.data[:,0])) 
print('x1?: ',list(raw1.x1()))    

print('tag?: ',list(raw1.data[:,7:])) 

print('proc_tag: ',list(raw1.proc_tag()))  
print('part_tag: ',list(raw1.part_tag()))  


# File 2
fname2 = sys.argv[2]
fpath2, fext2 = os.path.splitext(fname2)

print('Reading', fname2)
if any(fext2 == s for s in h5_ext_list):
	raw2 = H5RAW3D(fname2)
elif any(fext2 == s for s in ascii_ext_list):
	raw2 = CSVRAW3D(fname2)
else:
	print('ERROR:\t\tExtension of file:', str(fname2))
	print('\t\t\tis not supported!')
	print('Allowed hdf5-file extensions: ',list(h5_ext_list))
	print('Allowed csv-file extensions: ',list(ascii_ext_list))
	sys.exit()


raw2.printhdf5keys()

print('x1?: ',list(raw2.data[:,0])) 
print('x1?: ',list(raw2.x1()))    

print('tag?: ',list(raw2.data[:,7:])) 
print('tag?: ',list(raw2.tag()[:,1])) 

print('part_tag: ',list(raw2.part_tag()))    	    
