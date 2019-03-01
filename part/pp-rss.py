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
    parser.add_argument(  "--Nstride",
                          type=int,
                          action='store',
                          dest="Nstride",
                          metavar="NSTRIDE",
                          default=1,
                          help='Processing only every nth file.')    
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



def comp_slices(file, Nbins, zeta_range, cellvol, order, crossterms, showtimings):
    raw = HiRAW(file)
    raw.read_attrs()
    raw.read_data(verbose=False)

    time = raw.get_time()

    slices = pp_raw_ana.Slices(raw, nbins=Nbins, zeta_range=zeta_range, cellvol=cellvol)

    slices.calc_moments(order=order, crossterms=crossterms, showtimings=showtimings)

    return time, slices


def main():

    parser = ps_parseargs()

    args = parser.parse_args()

    mom_order = args.mom_order
    crossterms = args.crossterms

    h5fl = H5FList(args.path, h5ftype='raw')
    flist = h5fl.get(verbose=False, stride=args.Nstride)
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
        slices = pp_raw_ana.Slices(raw, zeta_range=zeta_range)          
        Nbins = slices.nbins
    else:
        Nbins = args.Nbins
        zeta_range = args.zeta_range    

    sys.stdout.write('There are %i raw files to process...\n' % Nfiles)
    sys.stdout.flush()

    sm = SliceMoms()
    sm.alloc(Nzeta = Nbins, Nt = Nfiles, order = mom_order, with_2nd_order_xterms = crossterms)

    Nworkers = 2

    print('Nfiles = %d' % Nfiles)

    for i in range(0, Nfiles, Nworkers):
        ltime = []
        lslices = []

        if (i + Nworkers < Nfiles):
            Nsubloops = Nworkers
        else:
            Nsubloops = Nfiles - i

        for j in range(0,Nsubloops):          
            idx = i+j

            print('idx = %d' % idx)

            file = flist[idx]

            sys.stdout.write('Processing: %s\t(%i/%i)\n' % (file, (idx+1), Nfiles))
            sys.stdout.flush()  

            time, slices = comp_slices( file=file, \
                                        Nbins=Nbins, \
                                        zeta_range=zeta_range, \
                                        cellvol=cellvol, \
                                        order=mom_order, \
                                        crossterms=crossterms, \
                                        showtimings=args.timings)
            ltime.append(time)
            lslices.append(slices)


        for j in range(0,Nsubloops):
            idx = i+j
            sm.set_time(time,idx)

            sm.set_at_nt(lslices[j].charge,idx)

            sm.set_at_nt(lslices[j].avgx1,idx,x1=1)
            sm.set_at_nt(lslices[j].avgx2,idx,x2=1)
            sm.set_at_nt(lslices[j].avgx3,idx,x3=1)

            sm.set_at_nt(lslices[j].avgp1,idx,p1=1)
            sm.set_at_nt(lslices[j].avgp2,idx,p2=1)
            sm.set_at_nt(lslices[j].avgp3,idx,p3=1)

            if mom_order>1:
                sm.set_at_nt(lslices[j].avgx1sq,idx,x1=2)
                sm.set_at_nt(lslices[j].avgx2sq,idx,x2=2)
                sm.set_at_nt(lslices[j].avgx3sq,idx,x3=2)
                
                sm.set_at_nt(lslices[j].avgp1sq,idx,p1=2)
                sm.set_at_nt(lslices[j].avgp2sq,idx,p2=2)
                sm.set_at_nt(lslices[j].avgp3sq,idx,p3=2)

                sm.set_at_nt(lslices[j].avgx1p1,idx,x1=1,p1=1)
                sm.set_at_nt(lslices[j].avgx2p2,idx,x2=1,p2=1)
                sm.set_at_nt(lslices[j].avgx3p3,idx,x3=1,p3=1)

                if crossterms:
                    sm.set_at_nt(lslices[j].avgx1x2,idx,x1=1,x2=1)
                    sm.set_at_nt(lslices[j].avgx1x3,idx,x1=1,x3=1)
                    sm.set_at_nt(lslices[j].avgx1p2,idx,x1=1,p2=1)
                    sm.set_at_nt(lslices[j].avgx1p3,idx,x1=1,p3=1)

                    sm.set_at_nt(lslices[j].avgx2x3,idx,x2=1,x3=1)
                    sm.set_at_nt(lslices[j].avgx2p1,idx,x2=1,p1=1)
                    sm.set_at_nt(lslices[j].avgx2p3,idx,x2=1,p3=1)

                    sm.set_at_nt(lslices[j].avgp1p2,idx,p1=1,p2=1)
                    sm.set_at_nt(lslices[j].avgp1p3,idx,p1=1,p3=1)
                    
                    sm.set_at_nt(lslices[j].avgp2p3,idx,p2=1,p3=1)

            if mom_order>2:
                sm.set_at_nt(lslices[j].avgx1cube,idx,x1=3)
                sm.set_at_nt(lslices[j].avgx2cube,idx,x2=3)
                sm.set_at_nt(lslices[j].avgx3cube,idx,x3=3)

                sm.set_at_nt(lslices[j].avgp1cube,idx,p1=3)
                sm.set_at_nt(lslices[j].avgp2cube,idx,p2=3)
                sm.set_at_nt(lslices[j].avgp3cube,idx,p3=3)

                sm.set_at_nt(lslices[j].avgx1sqp1,idx,x1=2,p1=1)
                sm.set_at_nt(lslices[j].avgx2sqp2,idx,x2=2,p2=1)
                sm.set_at_nt(lslices[j].avgx3sqp3,idx,x3=2,p3=1)

                sm.set_at_nt(lslices[j].avgx1p1sq,idx,x1=1,p1=2)
                sm.set_at_nt(lslices[j].avgx2p2sq,idx,x2=1,p2=2)
                sm.set_at_nt(lslices[j].avgx3p3sq,idx,x3=1,p3=2)

            if mom_order>3:
                sm.set_at_nt(lslices[j].avgx1quar,idx,x1=4)
                sm.set_at_nt(lslices[j].avgx2quar,idx,x2=4)
                sm.set_at_nt(lslices[j].avgx3quar,idx,x3=4)

                sm.set_at_nt(lslices[j].avgp1quar,idx,p1=4)
                sm.set_at_nt(lslices[j].avgp2quar,idx,p2=4)
                sm.set_at_nt(lslices[j].avgp3quar,idx,p3=4)

                sm.set_at_nt(lslices[j].avgx1sqp1sq,idx,x1=2,p1=2)
                sm.set_at_nt(lslices[j].avgx2sqp2sq,idx,x2=2,p2=2)
                sm.set_at_nt(lslices[j].avgx3sqp3sq,idx,x3=2,p3=2)

    sm.set_zeta_array(lslices[-1].centers)

    savepath = mkdirs_if_nexist(args.savepath)

    h5savepathname = savepath + '/' + args.save_name

    sys.stdout.write('Saving to file: %s\n' % (h5savepathname))
    sys.stdout.flush()

    sm.write(h5savepathname)

    sys.stdout.write('Done!\n')
    sys.stdout.flush()

if __name__ == "__main__":
    main()
