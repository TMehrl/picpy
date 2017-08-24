#!/usr/bin/env python3
# This script may be executed like this:
# nohup raw-slice-series.py <DATA> 1> rss.out 2> rss.err &

import numpy as np
import os
from argparse import ArgumentParser
#from argparse import OptionGroup
import math
import sys
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import picdefs
from h5dat import HiRAW
from h5dat import DIR
import ps_ana
import h5py

# Parse defaults/definitions
class parsedefaults:
  save_name = 'slice-avgs.h5'
  savepath = './'
  raw_ident_str = 'raw_beam_'
  mom_order = 2
  Nbins = 256
  Nfiles = 0
  crossterms = False

def ps_parseargs():

  usg = 'Usage: %prog [options] <file or path>'
  
  desc='This is the picpy postprocessing tool.'
  
  parser = ArgumentParser(usage=usg, description=desc)
  parser.add_argument(  'path', 
                        metavar='PATH',
                        help='Path to raw outputs.')
  parser.add_argument(  '-v', '--verbose',
                        action='store_true', 
                        dest='verbose', 
                        default=True,
                        help = 'Print info (Default).')
  parser.add_argument(  '-q', '--quiet',
                        action='store_false', 
                        dest='verbose',
                        help = 'Don''t print info.')
  parser.add_argument(  "-s", "--save-path", 
                        dest="savepath",
                        metavar="PATH",
                        default=parsedefaults.savepath,
                        help = 'Path to which generated files will be saved. '
                               '(Default: "%s")' % parsedefaults.savepath)
  parser.add_argument(  "--raw-istr", 
                        dest="raw_ident_str",
                        metavar="RAWIdentstr",
                        default=parsedefaults.raw_ident_str,
                        help = 'Identification string for beam raw file. '
                              '(Default: "%s")' % parsedefaults.raw_ident_str)                            
  parser.add_argument(  "-n", "--save-name", 
                        dest="save_name",
                        metavar="NAME",
                        default=parsedefaults.save_name,
                        help = 'Define customized output filename.')  
  parser.add_argument(  "-o", "--mom-order", 
                        type=int,                  
                        action='store',
                        dest="mom_order",
                        metavar="MOMORDER",
                        choices=[1, 2, 3,],
                        default=parsedefaults.mom_order,
                        help='Order of moment evaluation (Default: 2).'
                        ) 
#                        '(Default: %i).' % parsedefaults.mom_order)                                      
  parser.add_argument(  "--Nfiles",
                        type=int, 
                        action='store',
                        dest="Nfiles",
                        metavar="NFILES",
                        default=parsedefaults.Nfiles,
                        help='Number of files to analyze.')
  parser.add_argument(  "--Nbins",
                        type=int,
                        action='store',
                        dest="Nbins",
                        metavar="Nbins",
                        default=parsedefaults.Nbins,
                        help= 'Number of bins. (Default: %i)' % parsedefaults.Nbins)
  parser.add_argument(  "--xterms",
                        action="store_true", 
                        dest="crossterms", 
                        default=parsedefaults.crossterms,
                        help = 'Compute averages of crossterms, '
                        'e.g. <x1p3> (Default: %s).' % parsedefaults.crossterms)                        
#   parser.add_argument(  "-c", "--code", 
#                       type='choice',
#                       action='store',
#                       dest="piccode",
#                       metavar="CODE",
#                       choices = [picdefs.code.hipace, picdefs.code.osiris,],
#                       default = picdefs.code.hipace,
#                       help= "PIC code which was used to generate files (Default: '%s')." 
#                             % picdefs.code.hipace)
#   parser.add_argument(  "-d", "--dim", 
#                       type='choice',
#                       action='store',
#                       dest="dimensionality",
#                       metavar="DIM",
#                       choices=[1, 2, 3,],
#                       default=3,
#                       help= 'Dimensionality of PIC simulation
#                            (Default: 3).')                                             
#   group = OptionGroup(parser, "Options for beam-phase-space (RAW) files",
#                       "These are options for beam-phase-space (RAW) files")
#   group.add_argument("-g", action="store_true", help="Group option.")
#   parser.add_argument_group(group)

#   group = OptionGroup(parser, "Options for grid files",
#                       "These are options for grid files")
#   group.add_argument("-g", action="store_true", help="Group option.")
#   parser.add_argument_group(group)
  
  return parser



def main():
  
  parser = ps_parseargs()

  args = parser.parse_args()
  
  Nbins = args.Nbins
  mom_order = args.mom_order
  crossterms = args.crossterms

  dir = DIR(args.path)
  dir.list_files(args.raw_ident_str)

  if (dir.nf > 0):
    if (args.Nfiles == parsedefaults.Nfiles):
      Nfiles = dir.nf
    elif int(args.Nfiles) <= dir.nf:
      Nfiles = int(args.Nfiles)
    else:
      sys.stderr('Error: Nfiles cannot be greater than the actual number of files!')
  else :
    print(  'Error:\tNo phase space (raw) files in directory!\n' + 
            ('\tCheck the used path (currently: "%s")\n' % args.savepath) + 
            ('\tCheck the used RAW identification string (currently: "%s")' % args.raw_ident_str)
          )
    sys.exit()     

  sys.stdout.write('There are %i raw files to process...\n' % Nfiles)
  sys.stdout.flush()

  time_array = np.zeros(Nfiles, dtype=np.float32)
  avgx1 = np.zeros((Nfiles, Nbins), dtype=np.float32)
  avgx2 = np.zeros((Nfiles, Nbins), dtype=np.float32)
  avgx3 = np.zeros((Nfiles, Nbins), dtype=np.float32)    
  avgp1 = np.zeros((Nfiles, Nbins), dtype=np.float32)
  avgp2 = np.zeros((Nfiles, Nbins), dtype=np.float32)
  avgp3 = np.zeros((Nfiles, Nbins), dtype=np.float32)
  
  if mom_order>1:
    avgx1sq = np.zeros((Nfiles, Nbins), dtype=np.float32)
    avgx2sq = np.zeros((Nfiles, Nbins), dtype=np.float32)
    avgx3sq = np.zeros((Nfiles, Nbins), dtype=np.float32)    
    avgp1sq = np.zeros((Nfiles, Nbins), dtype=np.float32)
    avgp2sq = np.zeros((Nfiles, Nbins), dtype=np.float32)
    avgp3sq = np.zeros((Nfiles, Nbins), dtype=np.float32)  
    avgx1p1 = np.zeros((Nfiles, Nbins), dtype=np.float32)
    avgx2p2 = np.zeros((Nfiles, Nbins), dtype=np.float32)
    avgx3p3 = np.zeros((Nfiles, Nbins), dtype=np.float32)
    
    if crossterms:
      avgx1x2 = np.zeros((Nfiles, Nbins), dtype=np.float32)
      avgx3x1 = np.zeros((Nfiles, Nbins), dtype=np.float32)
      avgx2x3 = np.zeros((Nfiles, Nbins), dtype=np.float32)   
      avgp1p2 = np.zeros((Nfiles, Nbins), dtype=np.float32)
      avgp3p1 = np.zeros((Nfiles, Nbins), dtype=np.float32)
      avgp2p3 = np.zeros((Nfiles, Nbins), dtype=np.float32)             
      avgx1p2 = np.zeros((Nfiles, Nbins), dtype=np.float32)
      avgx1p3 = np.zeros((Nfiles, Nbins), dtype=np.float32)
      avgx2p1 = np.zeros((Nfiles, Nbins), dtype=np.float32) 
      avgx2p3 = np.zeros((Nfiles, Nbins), dtype=np.float32)       
      avgx3p1 = np.zeros((Nfiles, Nbins), dtype=np.float32)        
      avgx3p2 = np.zeros((Nfiles, Nbins), dtype=np.float32)
 
  if mom_order>2:
    avgx1cube = np.zeros((Nfiles, Nbins), dtype=np.float32)
    avgx2cube = np.zeros((Nfiles, Nbins), dtype=np.float32)
    avgx3cube = np.zeros((Nfiles, Nbins), dtype=np.float32)    
    avgp1cube = np.zeros((Nfiles, Nbins), dtype=np.float32)
    avgp2cube = np.zeros((Nfiles, Nbins), dtype=np.float32)
    avgp3cube = np.zeros((Nfiles, Nbins), dtype=np.float32)
    avgx1sqp1 = np.zeros((Nfiles, Nbins), dtype=np.float32)  
    avgx2sqp2 = np.zeros((Nfiles, Nbins), dtype=np.float32)
    avgx3sqp3 = np.zeros((Nfiles, Nbins), dtype=np.float32)
    avgx1p1sq = np.zeros((Nfiles, Nbins), dtype=np.float32)  
    avgx2p2sq = np.zeros((Nfiles, Nbins), dtype=np.float32)
    avgx3p3sq = np.zeros((Nfiles, Nbins), dtype=np.float32)

    
    if crossterms:
      print('Third order crossterms not yet implemented!') 

        
  for i in range(0, Nfiles):
    sys.stdout.write('Processing: %s\t(%i/%i)\n' % (dir.filepath(i), i+1, Nfiles))
    sys.stdout.flush()
    
    raw = HiRAW(dir.filepath(i))
    raw.read_attrs()
    raw.read_data()
    
    time_array[i] = raw.time
    slices = ps_ana.Slices(raw, nbins=Nbins)
    slices.calc_moments(order = mom_order, crossterms=crossterms)
    avgx1[i,:] = slices.avgx1      
    avgx2[i,:] = slices.avgx2
    avgx3[i,:] = slices.avgx3
    avgp1[i,:] = slices.avgp1      
    avgp2[i,:] = slices.avgp2
    avgp3[i,:] = slices.avgp3  
    
    if mom_order>1:
      avgx1sq[i,:] = slices.avgx1sq      
      avgx2sq[i,:] = slices.avgx2sq
      avgx3sq[i,:] = slices.avgx3sq
      avgp1sq[i,:] = slices.avgp1sq    
      avgp2sq[i,:] = slices.avgp2sq
      avgp3sq[i,:] = slices.avgp3sq
      avgx1p1[i,:] = slices.avgx1p1    
      avgx2p2[i,:] = slices.avgx2p2
      avgx3p3[i,:] = slices.avgx3p3
      
      if crossterms:
        avgx1x2[i,:] = slices.avgx1x2 
        avgx3x1[i,:] = slices.avgx3x1 
        avgx2x3[i,:] = slices.avgx2x3    
        avgp1p2[i,:] = slices.avgp1p2 
        avgp3p1[i,:] = slices.avgp3p1 
        avgp2p3[i,:] = slices.avgp2p3              
        avgx1p2[i,:] = slices.avgx1p2 
        avgx1p3[i,:] = slices.avgx1p3 
        avgx2p1[i,:] = slices.avgx2p1  
        avgx2p3[i,:] = slices.avgx2p3        
        avgx3p1[i,:] = slices.avgx3p1         
        avgx3p2[i,:] = slices.avgx3p2       
    
    if mom_order>2:
      avgx1cube[i,:] = slices.avgx1cube      
      avgx2cube[i,:] = slices.avgx2cube
      avgx3cube[i,:] = slices.avgx3cube
      avgp1cube[i,:] = slices.avgp1cube    
      avgp2cube[i,:] = slices.avgp2cube
      avgp3cube[i,:] = slices.avgp3cube
      avgx1sqp1[i,:] = slices.avgx1sqp1    
      avgx2sqp2[i,:] = slices.avgx2sqp2
      avgx3sqp3[i,:] = slices.avgx3sqp3
      avgx1p1sq[i,:] = slices.avgx1p1sq    
      avgx2p2sq[i,:] = slices.avgx2p2sq
      avgx3p3sq[i,:] = slices.avgx3p3sq
            
    
  zeta_array = slices.centers
  
  h5savepathname = args.savepath + '/' + args.save_name
  
  sys.stdout.write('Saving to file: %s\n' % (h5savepathname))
  sys.stdout.flush()

  h5f = h5py.File(h5savepathname, "w")
  dset_zeta_array = h5f.create_dataset( "zeta_array", data = zeta_array )
  dset_time_array = h5f.create_dataset( "time_array", data = time_array )
  dset_avgx1 = h5f.create_dataset(  "avgx1", data = avgx1 )
  dset_avgx2 = h5f.create_dataset(  "avgx2", data = avgx2 )
  dset_avgx3 = h5f.create_dataset(  "avgx3", data = avgx3 )
  dset_avgp1 = h5f.create_dataset(  "avgp1", data = avgp1 )
  dset_avgp2 = h5f.create_dataset(  "avgp2", data = avgp2 )
  dset_avgp3 = h5f.create_dataset(  "avgp3", data = avgp3 )

  if mom_order>1:
    dset_avgx1sq = h5f.create_dataset(  "avgx1sq", data = avgx1sq )
    dset_avgx2sq = h5f.create_dataset(  "avgx2sq", data = avgx2sq )
    dset_avgx3sq = h5f.create_dataset(  "avgx3sq", data = avgx3sq )
    dset_avgp1sq = h5f.create_dataset(  "avgp1sq", data = avgp1sq )
    dset_avgp2sq = h5f.create_dataset(  "avgp2sq", data = avgp2sq )
    dset_avgp3sq = h5f.create_dataset(  "avgp3sq", data = avgp3sq )
    dset_avgx1p1 = h5f.create_dataset(  "avgx1p1", data = avgx1p1 )
    dset_avgx2p2 = h5f.create_dataset(  "avgx2p2", data = avgx2p2 )
    dset_avgx3p3 = h5f.create_dataset(  "avgx3p3", data = avgx3p3 )

    if crossterms:
      dset_avgx1x2 = h5f.create_dataset(  "avgx1x2", data = avgx1x2 ) 
      dset_avgx3x1 = h5f.create_dataset(  "avgx3x1", data = avgx3x1 )
      dset_avgx2x3 = h5f.create_dataset(  "avgx2x3", data = avgx2x3 )   
      dset_avgp1p2 = h5f.create_dataset(  "avgp1p2", data = avgp1p2 )
      dset_avgp3p1 = h5f.create_dataset(  "avgp3p1", data = avgp3p1 ) 
      dset_avgp2p3 = h5f.create_dataset(  "avgp2p3", data = avgp2p3 )              
      dset_avgx1p2 = h5f.create_dataset(  "avgx1p2", data = avgx1p2 )
      dset_avgx1p3 = h5f.create_dataset(  "avgx1p3", data = avgx1p3 )
      dset_avgx2p1 = h5f.create_dataset(  "avgx2p1", data = avgx2p1 )
      dset_avgx2p3 = h5f.create_dataset(  "avgx2p3", data = avgx2p3 )     
      dset_avgx3p1 = h5f.create_dataset(  "avgx3p1", data = avgx3p1 )      
      dset_avgx3p2 = h5f.create_dataset(  "avgx3p2", data = avgx3p2 )    

  if mom_order>2:
    dset_avgx1cube = h5f.create_dataset(  "avgx1cube", data = avgx1cube )
    dset_avgx2cube = h5f.create_dataset(  "avgx2cube", data = avgx2cube )
    dset_avgx3cube = h5f.create_dataset(  "avgx3cube", data = avgx3cube )
    dset_avgp1cube = h5f.create_dataset(  "avgp1cube", data = avgp1cube ) 
    dset_avgp2cube = h5f.create_dataset(  "avgp2cube", data = avgp2cube )
    dset_avgp3cube = h5f.create_dataset(  "avgp3cube", data = avgp3cube )
    dset_avgx1sqp1 = h5f.create_dataset(  "avgx1sqp1", data = avgx1sqp1 )  
    dset_avgx2sqp2 = h5f.create_dataset(  "avgx2sqp2", data = avgx2sqp2 )
    dset_avgx3sqp3 = h5f.create_dataset(  "avgx3sqp3", data = avgx3sqp3 )
    dset_avgx1p1sq = h5f.create_dataset(  "avgx1p1sq", data = avgx1p1sq ) 
    dset_avgx2p2sq = h5f.create_dataset(  "avgx2p2sq", data = avgx2p2sq )
    dset_avgx3p3sq = h5f.create_dataset(  "avgx3p3sq", data = avgx3p3sq )
      
  h5f.close() 
                                              
  sys.stdout.write('Done!\n')
  sys.stdout.flush()
  
if __name__ == "__main__":
    main()