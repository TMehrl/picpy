#!/usr/bin/env python3

import argparse
import plot_g3d

def main():
  parser = plot_g3d.parser( ptype='slice' )
  args = parser.parse_args()
  plot_g3d.plotfiles( args, ptype='slice' )
   
if __name__ == "__main__":
    main()