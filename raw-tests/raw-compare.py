#!/usr/bin/env python3.5

import h5py
import numpy as np
import sys
from optparse import OptionParser
from scipy import constants
import os
import csv
import math

# Global
h5_ext_list = ['.h5','.hdf5']
ascii_ext_list = ['.csv','.txt','.ascii']

def read_file(fname):

	fpath, fext = os.path.splitext(fname)

	print('Reading', fname)
	if any(fext == s for s in h5_ext_list) or any(fext == s for s in ascii_ext_list):
		raw = RAW3D(fname)
	else:
		print('ERROR:\tExtension of file "%s" is not supported!' %fname)
		print('\tAllowed hdf5-file extensions: ',list(h5_ext_list))
		print('\tAllowed csv-file extensions: ',list(ascii_ext_list))
		sys.exit()
	print('Reading complete.')
	return raw


class RAW3D:
  def __init__(self, file):
  	fpath, fext = os.path.splitext(file)
  	self.filename = file
  	if any(fext == s for s in h5_ext_list):
  		self.dattype = 'hdf5'
  		with h5py.File(file,'r') as hf:
  			self.hdf5keys = list(hf.keys())
  			self.x1 = np.array(hf.get('x1'))
  			self.x2 = np.array(hf.get('x2'))
  			self.x3 = np.array(hf.get('x3'))
  			self.p1 = np.array(hf.get('p1'))
  			self.p2 = np.array(hf.get('p2'))
  			self.p3 = np.array(hf.get('p3'))
  			self.q  = np.array(hf.get('q'))
  			self.tag  = np.array(hf.get('tag'))
  			hf.close()
  	elif any(fext == s for s in ascii_ext_list):
  		self.dattype = 'csv'
  		with open(file, 'rt') as cf:
  			csvstrdata = csv.reader(cf, delimiter='\t', quoting=csv.QUOTE_NONNUMERIC)
  			self.data = np.array(list(csvstrdata)).astype('float')
  			cf.close()
  		self.x1 = self.data[:,0]
  		self.x2 = self.data[:,1]
  		self.x3 = self.data[:,2]
  		self.p1 = self.data[:,3]
  		self.p2 = self.data[:,4]
  		self.p3 = self.data[:,5]
  		self.q = self.data[:,6]
  		self.tag  = np.zeros((len(self.data[:,7]),2))
  		self.tag[:,0]  = np.ones(len(self.data[:,7]))
  		self.tag[:,1] = np.array(self.data[:,7]).astype(int)
  					
  def proc_tag(self):
  	return self.tag[:,0]  	
  def part_tag(self):
  	return self.tag[:,1]			  
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
  def conv_to_si(self, n0):
  	k_p = math.sqrt(n0 * constants.e**2/(constants.epsilon_0 * constants.m_e))\
  	/constants.c
  	self.x1 = self.x1*k_p
  	self.x2 = self.x2*k_p
  	self.x3 = self.x3*k_p

def main():
	NUM_ARGS = 2
	PLASMA_DEN_DEFAULT_STR = "0.0"

	usage = "usage: %prog [options] <beam-file> <raw-file>"
	parser = OptionParser(usage=usage)
	parser.add_option("-n", "--plasma-density", dest="n0",
					  metavar="DENSITY", default=PLASMA_DEN_DEFAULT_STR,
					  help="Plasma density n0 in m^-3. "
					  "If set, original beam-file is assumed in Si units. "
					  "Write exponentiation as 'XeY'!")

	(options, args) = parser.parse_args()

	if len(args) != 2:
		parser.error("This script requires two arguments!")	

	# File 1
	raw1=read_file(args[0])
	# File 2
	raw2=read_file(args[1])

	if not options.n0 == PLASMA_DEN_DEFAULT_STR:
		n0 = float(options.n0)
		print('raw2.x1=',raw2.x1,)
		raw2.conv_to_si(n0)
		print('raw2.x1=',raw2.x1,)

	raw1.sort()
	raw2.sort()

	if (max(abs(np.subtract(raw1.part_tag(), raw2.part_tag())))==0):
		print( 'Tags correspond.')
	else:
		print( 'Tags do NOT correspond!')

	print( 'Max. x1-diff: %.15f' % max(abs(np.subtract(raw1.x1, raw2.x1))) )
	print( 'Max. x2-diff: %.15f' % max(abs(np.subtract(raw1.x2, raw2.x2))) )
	print( 'Max. x3-diff: %.15f' % max(abs(np.subtract(raw1.x3, raw2.x3))) )
	print( 'Max. p1-diff: %.15f' % max(abs(np.subtract(raw1.p1, raw2.p1))) )
	print( 'Max. p2-diff: %.15f' % max(abs(np.subtract(raw1.p2, raw2.p2))) )
	print( 'Max. p3-diff: %.15f' % max(abs(np.subtract(raw1.p3, raw2.p3))) )
	print( 'Max. q-diff: %.15f' % max(abs(np.subtract(raw1.q, raw2.q))) )


if __name__ == "__main__":
    main()