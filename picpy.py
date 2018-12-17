#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
import argcomplete

# get picpy package root path
pp_path = os.path.dirname(os.path.abspath( __file__ ))

# set absolute paths of directories 
# with the different types of scripts
part_path = pp_path + '/part/'
grid_path = pp_path + '/grid/'
special_path = pp_path + '/special/'

# A dictionary of all available scripts/tools
# and their respective absolute paths
# (This may at some point also automatically
# be generated by browsing the directories)
scripts = {
  "g3d": grid_path + "pp-g3d.py",
  "g3dline" : grid_path + "pp-g3dline.py",
  "rss": part_path + "pp-rss.py",
  "rss-plot": part_path + "pp-rss-plot.py",
  "raw2hist": part_path + "pp-raw2hist.py",
  "3trackingbin" : special_path + "pp-3trackingbin.py",
  "binLine" : special_path + "pp-binLine.py",
  "binSlab" : special_path + "pp-binSlab.py",
  "blowout-geo" : special_path + "pp-blowout-geo.py",
  "calc_Ez_slope" : special_path + "pp-calc_Ez_slope.py",
  "comparison" : special_path + "pp-comparison.py",
  "equilib" : special_path + "pp-equilib.py",
  "ex_ioniz_dens_lineout" : special_path + "pp-ex_ioniz_dens_lineout.py",
  "gauss_fit" : special_path + "pp-gauss_fit.py",
  "h5plot" : special_path + "pp-h5plot.py",
  "ionization" : special_path + "pp-ionization.py",
  "lvt" : special_path + "pp-lvt.py",
  "mkmov" : special_path + "pp-mkmov.py",
  "raw-test" : special_path + "pp-raw-test.py",
  "tracking" : special_path + "pp-tracking.py",
  "tracking3d" : special_path + "pp-tracking3d.py",
  "trackingbin" : special_path + "pp-trackingbin.py"
}

def main(**args):
    """
    picpy main function.

    Function calls scripts according to args.

    Parameters
    ----------
    args : dict
        args contains the information which script
        (keyword 'script') is to be called with what
        arguments (keyword 'args'). 
        If args['bkg'] == True, the script is called
        in background with the POSIX command 'nohup'. 
    """  
    exelist = [scripts[args['script']]] + args['args']
    if args['bkg']:
        print('Executing %s in background.' \
                '\nWriting stdout to "pp.out" and stderr to "pp.err"'\
                % args['script'])
        subprocess.Popen(['nohup'] + exelist,
                         stdout=open('pp.out', 'w'),
                         stderr=open('pp.err', 'a'))
    else:
        subprocess.call(exelist)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('script', 
                        choices=list(scripts.keys()))
    parser.add_argument('-b', '--bkg',
                        dest="bkg",
                        action="store_true",
                        default=False,
                        help="Run script in background" 
                            "(Default: %(default)s).")
    parser.add_argument('args', 
                        nargs=argparse.REMAINDER )
    argcomplete.autocomplete(parser)
    args = parser.parse_args()
    main(**vars(args))
