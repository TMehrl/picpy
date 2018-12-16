#!/usr/bin/env python3

import os
import sys
import shutil
import subprocess
import argparse


def ps_parseopts():

    desc="""Movie generation tool."""

    parser = argparse.ArgumentParser(description=desc)
    #parser = OptionParser(usage=usg, description=desc)
    parser.add_argument(  'path',
                          metavar = 'PATH',
                          nargs = '*',                          
                          help = 'Path to frames.')    
    parser.add_argument(  "-v", "--verbose",
                        action="store_true",
                        dest="verbose",
                        default=True,
                        help = "Print info (Default: %(default)s).")
    parser.add_argument(  "-q", "--quiet",
                        action="store_false",
                        dest="verbose",
                        help = "Don't print info.")
    parser.add_argument(  "--inf",
                          action='store',
                          dest="input_format",
                          metavar="INPUT_FORMAT",
                          choices=[ 'png',
                                    'jpeg',
                                    'eps',],
                          default='png',
                          help= """Format of input frames (Default:  %(default)s).""")
    parser.add_argument(  "--movname",
                          action='store',
                          dest="movname",
                          metavar="MOVIE_NAME",
                          default='mov',
                          help= """Name of output movie (Default:  %(default)s).""") 
    parser.add_argument(  "--framerate",
                          action='store',
                          dest="framerate",
                          metavar="FRAMERATE",
                          type=int,                          
                          default=10,
                          help= """Framerate (Default:  %(default)s).""")                                                            
    parser.add_argument(  "-s", "--save-path",
                          action="store",
                          dest="savepath",
                          metavar="PATH",
                          default='.',
                          help = """Path to which generated movie will be saved (Default: %(default)s)""") 
    return parser  


def get_filelist(path):
    
    if len(path) == 0:
        print('ERROR: Expecting path argument!')
        sys.exit()        
    elif len(path) == 1:
        if os.path.isdir(path[0]):
            print('Collecting all files in %s' % path[0])
            files = []
            for (dirpath, dirnames, filenames) in os.walk(path[0]):
                for file in filenames:
                    files.append(os.path.join(dirpath, file))
        else:
            print('ERROR: Expecting directory for single path argument!')
            sys.exit() 
    else:
        files = path    

    return files


def fametime(file):
    filename, file_extension = os.path.splitext(file)
    return(filename[-6:])

def collect_sort_frames(filelist, ext = 'png'):

    framelist = []

    for file in filelist:
        filename, file_extension = os.path.splitext(file)    
        if file_extension == ('.' + ext):
            framelist.append(file)
        else:
            print('Skipping %s! Wrong file format.' % file)
            print('Expecting %s-files!' % ext)

    return sorted(framelist, key = fametime)

def main():

    parser = ps_parseopts()

    args = parser.parse_args()  

    temp_frame_path="frames"
    frame_prefix = "frame"

    filelist = get_filelist(args.path)
    framelist = collect_sort_frames(filelist)

    os.mkdir('./' + temp_frame_path)

    for i in range(len(framelist)):
        frame = framelist[i]
        shutil.copy(frame, '%s/%s_%06i.%s' % (temp_frame_path, frame_prefix, i, args.input_format))

    # Generating mp4 movie using ffmpeg
    makemov = subprocess.call(["ffmpeg", \
    "-framerate", ("%i" % args.framerate), \
    "-i", (temp_frame_path + "/" + frame_prefix + "_%06d." + args.input_format), \
    "-c:v", "libx264", \
    "-pix_fmt", "yuv420p", \
    args.savepath + "/" + args.movname + ".mp4"]) 

    shutil.rmtree(temp_frame_path)

if __name__ == '__main__':
    main()