#!/usr/bin/env python3

import os
import sys
import argparse
import argcomplete
import subprocess

pp_path = os.path.dirname(os.path.abspath( __file__ ))

part_path = pp_path + '/part/'
grid_path = pp_path + '/grid/'
special_path = pp_path + '/special/'

scripts =   {
  "g3d": (grid_path + "pp-g3d.py"),
  "rss": (part_path + "pp-rss.py"),
  "rss-plot": (part_path + "pp-rss-plot.py"),
  "raw2hist": (part_path + "pp-raw2hist.py")
}

def main(**args):
    exelist = [scripts[args['script']]] + args['args']
    print(scripts[args['script']]) 
    print(args['args'])
    print(exelist)
    subprocess.call(exelist)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('script', choices=list(scripts.keys()))
    parser.add_argument('args', nargs=argparse.REMAINDER )
    argcomplete.autocomplete(parser)
    args = parser.parse_args()
    main(**vars(args))