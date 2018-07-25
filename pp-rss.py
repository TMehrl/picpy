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
import pp_defs
from pp_h5dat import HiRAW
from pp_h5dat import H5FList
from pp_h5dat import SliceMoms
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
    parser.add_argument(  "-s", "--save-path",
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
        sys.exit()

    if args.Nfiles == None:
        Nfiles = len(flist)
    else:
        Nfiles = args.Nfiles

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
    sm.alloc(Nzeta = Nbins, Nt = Ntimesteps, order = mom_order)

    for i in range(0, Ntimesteps):
        file = flist[i * args.Nskip]
        sys.stdout.write('Processing: %s\t(%i/%i)\n' % (file, (i+1)*args.Nskip, Nfiles))
        sys.stdout.flush()

        raw = HiRAW(file)
        raw.read_attrs()
        raw.read_data()

        sm.time_array[i] = raw.get_time()
        slices = pp_raw_ana.Slices(raw, nbins=Nbins, zrange=zeta_range, cellvol=cellvol)

        slices.calc_moments(order = mom_order, crossterms=crossterms, timings=args.timings )
        sm.charge[i,:] = slices.charge
        sm.avgx1[i,:] = slices.avgx1
        sm.avgx2[i,:] = slices.avgx2
        sm.avgx3[i,:] = slices.avgx3
        sm.avgp1[i,:] = slices.avgp1
        sm.avgp2[i,:] = slices.avgp2
        sm.avgp3[i,:] = slices.avgp3

        if mom_order>1:
            sm.avgx1sq[i,:] = slices.avgx1sq
            sm.avgx2sq[i,:] = slices.avgx2sq
            sm.avgx3sq[i,:] = slices.avgx3sq
            sm.avgp1sq[i,:] = slices.avgp1sq
            sm.avgp2sq[i,:] = slices.avgp2sq
            sm.avgp3sq[i,:] = slices.avgp3sq
            sm.avgx1p1[i,:] = slices.avgx1p1
            sm.avgx2p2[i,:] = slices.avgx2p2
            sm.avgx3p3[i,:] = slices.avgx3p3

            if crossterms:
                sm.avgx1x2[i,:] = slices.avgx1x2
                sm.avgx3x1[i,:] = slices.avgx3x1
                sm.avgx2x3[i,:] = slices.avgx2x3
                sm.avgp1p2[i,:] = slices.avgp1p2
                sm.avgp3p1[i,:] = slices.avgp3p1
                sm.avgp2p3[i,:] = slices.avgp2p3
                sm.avgx1p2[i,:] = slices.avgx1p2
                sm.avgx1p3[i,:] = slices.avgx1p3
                sm.avgx2p1[i,:] = slices.avgx2p1
                sm.avgx2p3[i,:] = slices.avgx2p3
                sm.avgx3p1[i,:] = slices.avgx3p1
                sm.avgx3p2[i,:] = slices.avgx3p2

        if mom_order>2:
            sm.avgx1cube[i,:] = slices.avgx1cube
            sm.avgx2cube[i,:] = slices.avgx2cube
            sm.avgx3cube[i,:] = slices.avgx3cube
            sm.avgp1cube[i,:] = slices.avgp1cube
            sm.avgp2cube[i,:] = slices.avgp2cube
            sm.avgp3cube[i,:] = slices.avgp3cube
            sm.avgx1sqp1[i,:] = slices.avgx1sqp1
            sm.avgx2sqp2[i,:] = slices.avgx2sqp2
            sm.avgx3sqp3[i,:] = slices.avgx3sqp3
            sm.avgx1p1sq[i,:] = slices.avgx1p1sq
            sm.avgx2p2sq[i,:] = slices.avgx2p2sq
            sm.avgx3p3sq[i,:] = slices.avgx3p3sq


        if mom_order>3:
            sm.avgx1quar[i,:] = slices.avgx1quar
            sm.avgx2quar[i,:] = slices.avgx2quar
            sm.avgx3quar[i,:] = slices.avgx3quar
            sm.avgp1quar[i,:] = slices.avgp1quar
            sm.avgp2quar[i,:] = slices.avgp2quar
            sm.avgp3quar[i,:] = slices.avgp3quar
            sm.avgx1sqp1sq[i,:] = slices.avgx1sqp1sq
            sm.avgx2sqp2sq[i,:] = slices.avgx2sqp2sq
            sm.avgx3sqp3sq[i,:] = slices.avgx3sqp3sq

    sm.zeta_array = slices.centers

    h5savepathname = args.savepath + '/' + args.save_name

    sys.stdout.write('Saving to file: %s\n' % (h5savepathname))
    sys.stdout.flush()

    sm.write(h5savepathname)

    sys.stdout.write('Done!\n')
    sys.stdout.flush()

if __name__ == "__main__":
    main()
