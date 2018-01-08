#!/usr/bin/env python3

import plot_g3d

def main():
    parser = plot_g3d.parser( ptype='slice' )
    args = parser.parse_args()
   
    flist = plot_g3d.gen_filelist(args)

    for file in flist:
        g3d_p = plot_g3d.G3d_plot_slice(file, args)
        g3d_p.plot()

if __name__ == "__main__":
    main()