#!/usr/bin/env python3
# This script may be executed like this:
# nohup pp-rss.py DATA 1> rss.out 2> rss.err &

import os
import sys
import math
import argparse
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm

mypath = os.path.dirname(os.path.abspath( __file__ ))
incpath = os.path.split(mypath)[0] + '/inc'
sys.path.append(incpath)
import pp_defs
from pp_h5dat import HiRAW
from pp_h5dat import H5FList
from pp_h5dat import SliceMoms
from pp_h5dat import mkdirs_if_nexist
import pp_raw_ana

# Parse defaults/definitions
class parsedefaults:
    save_name = 'slice-avgs.h5'
    savepath = './'
    raw_ident_str = 'raw_beam_'
    mom_order = 2
    crossterms = False

def ps_parseargs():

    desc='This is the picpy postprocessing tool.'

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(  'path',
                          metavar = 'PATH',
                          nargs = '*',
                          help = 'Path to raw files.')
    parser.add_argument(  '-v', '--verbose',
                          action='store_true',
                          dest='verbose',
                          default=True,
                          help = 'Print info (Default).')
    parser.add_argument(  '-q', '--quiet',
                          action='store_false',
                          dest='verbose',
                          help = 'Don''t print info.')
    parser.add_argument(  "--save-path",
                          dest="savepath",
                          metavar="PATH",
                          default=parsedefaults.savepath,
                          help = 'Path to which generated files will be saved. '
                                 '(Default: "%(default)s")')
    parser.add_argument(  "--raw-istr",
                          dest="raw_ident_str",
                          metavar="RAWIdentstr",
                          default=parsedefaults.raw_ident_str,
                          help = 'Identification string for beam raw file. '
                                '(Default: "%(default)s")')
    parser.add_argument(  "--save-name",
                          dest="save_name",
                          metavar="NAME",
                          default=parsedefaults.save_name,
                          help = 'Define customized output filename.')
    parser.add_argument(  "-o", "--mom-order",
                          type=int,
                          action='store',
                          dest="mom_order",
                          metavar="MOMORDER",
                          choices=[1, 2, 3, 4,],
                          default=parsedefaults.mom_order,
                          help='Order of moment evaluation (Default: %(default)s).')
#                        '(Default: %i).' % parsedefaults.mom_order)
    parser.add_argument(  "--Nfiles",
                          type=int,
                          action='store',
                          dest="Nfiles",
                          metavar="NFILES",
                          default=None,
                          help='Number of files to analyze.')
    parser.add_argument(  "--Nskip",
                          type=int,
                          action='store',
                          dest="Nskip",
                          metavar="NSKIP",
                          default=1,
                          help='Skipping every nth file.')    
    parser.add_argument(  "--Nbins",
                          type=int,
                          action='store',
                          dest="Nbins",
                          metavar="Nbins",
                          default=None,
                          help= 'Number of bins.')
    parser.add_argument(  "--xterms",
                          action="store_true",
                          dest="crossterms",
                          default=parsedefaults.crossterms,
                          help = 'Compute averages of crossterms, '
                          'e.g. <x1p3> (Default: %s).' % parsedefaults.crossterms)
    parser.add_argument(  '--zeta-range',
                          help='zeta range',
                          action='store',
                          dest="zeta_range",
                          metavar=('ZETA_MIN', 'ZETA_MAX'),
                          nargs=2,
                          type=float,
                          default=None)
    parser.add_argument(  '-t', '--timings',
                          action='store_true',
                          dest='timings',
                          default=False,
                          help = 'Generate and print timings (default: %(default)s).')

#   parser.add_argument(  "-c", "--code",
#                       type='choice',
#                       action='store',
#                       dest="piccode",
#                       metavar="CODE",
#                       choices = [pp_defs.code.hipace, pp_defs.code.osiris,],
#                       default = pp_defs.code.hipace,
#                       help= "PIC code which was used to generate files (Default: '%s')."
#                             % pp_defs.code.hipace)
#   parser.add_argument(  "-d", "--dim",
#                       type='choice',
#                       action='store',
#                       dest="dimensionality",
#                       metavar="DIM",
#                       choices=[1, 2, 3,],
#                       default=3,
#                       help= 'Dimensionality of PIC simulation
#                            (Default: 3).')

    return parser



def main():

    parser = ps_parseargs()

    args = parser.parse_args()

    mom_order = args.mom_order
    crossterms = args.crossterms

    h5fl = H5FList(args.path, h5ftype='raw')
    flist = h5fl.get(verbose=False)
    if len(h5fl.get_uniques()) > 1:
        print('ERROR: Processing of multiple beams is not implemented yet!')
        print(h5fl.split_by_uniques())
        sys.exit(1)

    if args.Nfiles == None:
        Nfiles = len(flist)
    else:
        Nfiles = args.Nfiles

    if Nfiles < 1:
        print('No raw files selected!')
        print('Exiting...')
        sys.exit(1)        

    # Getting file information
    raw = HiRAW(flist[0])
    raw.read_attrs()
    cellvol = raw.get_dx(0) * raw.get_dx(1) * raw.get_dx(2)

    if args.Nbins == None and args.zeta_range == None:
        zeta_range = None
        slices = pp_raw_ana.Slices(raw)          
        Nbins = slices.nbins
    elif args.Nbins != None and args.zeta_range == None:
        Nbins = args.Nbins
        zeta_range = None
    elif args.Nbins == None and args.zeta_range != None:
        zeta_range = args.zeta_range
        slices = pp_raw_ana.Slices(raw, zrange=zeta_range)          
        Nbins = slices.nbins
    else:
        Nbins = args.Nbins
        zeta_range = args.zeta_range    

    sys.stdout.write('There are %i raw files to process...\n' % Nfiles)
    sys.stdout.flush()

    Ntimesteps = int( math.ceil(Nfiles/args.Nskip) )

    sm = SliceMoms()
    sm.alloc(Nzeta = Nbins, Nt = Ntimesteps, order = mom_order, with_2nd_order_xterms = crossterms)

    for i in range(0, Ntimesteps):
        file = flist[i * args.Nskip]
        sys.stdout.write('Processing: %s\t(%i/%i)\n' % (file, (i+1)*args.Nskip, Nfiles))
        sys.stdout.flush()

        raw = HiRAW(file)
        raw.read_attrs()
        raw.read_data(verbose=False)

        sm.set_time(raw.get_time(),i)
        slices = pp_raw_ana.Slices(raw, nbins=Nbins, zrange=zeta_range, cellvol=cellvol)

        slices.calc_moments(order = mom_order, crossterms=crossterms, showtimings=args.timings )

        sm.set_at_nt(slices.charge,i)

        sm.set_at_nt(slices.avgx1,i,x1=1)
        sm.set_at_nt(slices.avgx2,i,x2=1)
        sm.set_at_nt(slices.avgx3,i,x3=1)

        sm.set_at_nt(slices.avgp1,i,p1=1)
        sm.set_at_nt(slices.avgp2,i,p2=1)
        sm.set_at_nt(slices.avgp3,i,p3=1)

        if mom_order>1:
            sm.set_at_nt(slices.avgx1sq,i,x1=2)
            sm.set_at_nt(slices.avgx2sq,i,x2=2)
            sm.set_at_nt(slices.avgx3sq,i,x3=2)
            
            sm.set_at_nt(slices.avgp1sq,i,p1=2)
            sm.set_at_nt(slices.avgp2sq,i,p2=2)
            sm.set_at_nt(slices.avgp3sq,i,p3=2)

            sm.set_at_nt(slices.avgx1p1,i,x1=1,p1=1)
            sm.set_at_nt(slices.avgx2p2,i,x2=1,p2=1)
            sm.set_at_nt(slices.avgx2p2,i,x3=1,p3=1)

            if crossterms:
                sm.set_at_nt(slices.avgx1x2,i,x1=1,x2=1)
                sm.set_at_nt(slices.avgx3x1,i,x1=1,x3=1)
                sm.set_at_nt(slices.avgx3x1,i,x1=1,p2=1)
                sm.set_at_nt(slices.avgx3x1,i,x1=1,p3=1)

                sm.set_at_nt(slices.avgx2x3,i,x2=1,x3=1)
                sm.set_at_nt(slices.avgx2x3,i,x2=1,p1=1)
                sm.set_at_nt(slices.avgx2x3,i,x2=1,p3=1)

                sm.set_at_nt(slices.avgp1p2,i,p1=1,p2=1)
                sm.set_at_nt(slices.avgp3p1,i,p1=1,p3=1)
                
                sm.set_at_nt(slices.avgp2p3,i,p2=1,p3=1)


        if mom_order>2:
            sm.set_at_nt(slices.avgx1cube,i,x1=3)
            sm.set_at_nt(slices.avgx2cube,i,x2=3)
            sm.set_at_nt(slices.avgx3cube,i,x3=3)

            sm.set_at_nt(slices.avgp1cube,i,p1=3)
            sm.set_at_nt(slices.avgp2cube,i,p2=3)
            sm.set_at_nt(slices.avgp3cube,i,p3=3)

            sm.set_at_nt(slices.avgx1sqp1,i,x1=2,p1=1)
            sm.set_at_nt(slices.avgx2sqp2,i,x2=2,p2=1)
            sm.set_at_nt(slices.avgx3sqp3,i,x3=2,p3=1)

            sm.set_at_nt(slices.avgx1p1sq,i,x1=1,p1=2)
            sm.set_at_nt(slices.avgx2p2sq,i,x2=1,p2=2)
            sm.set_at_nt(slices.avgx3p3sq,i,x3=1,p3=2)

        if mom_order>3:
            sm.set_at_nt(slices.avgx1quar,i,x1=4)
            sm.set_at_nt(slices.avgx2quar,i,x2=4)
            sm.set_at_nt(slices.avgx3quar,i,x3=4)

            sm.set_at_nt(slices.avgp1quar,i,p1=4)
            sm.set_at_nt(slices.avgp2quar,i,p2=4)
            sm.set_at_nt(slices.avgp3quar,i,p3=4)

            sm.set_at_nt(slices.avgx1sqp1sq,i,x1=2,p1=2)
            sm.set_at_nt(slices.avgx2sqp2sq,i,x2=2,p2=2)
            sm.set_at_nt(slices.avgx3sqp3sq,i,x3=2,p3=2)

    sm.set_zeta_array(slices.centers)

    savepath = mkdirs_if_nexist(args.savepath)

    h5savepathname = savepath + '/' + args.save_name

    sys.stdout.write('Saving to file: %s\n' % (h5savepathname))
    sys.stdout.flush()

    sm.write(h5savepathname)

    sys.stdout.write('Done!\n')
    sys.stdout.flush()

if __name__ == "__main__":
    main()
