#!/usr/bin/env python3

import os
import sys
import importlib
import subprocess

argcomplete_spec = importlib.util.find_spec("argcomplete")
found_argcomplete = argcomplete_spec is not None
if found_argcomplete:
    import argcomplete

def main():
    """
    setup main function.

    Set environment variables for picpy
    execution and for argcomplete.
    """ 
    pp_path = os.path.dirname(os.path.abspath( __file__ ))

    if "PP_PATH" not in os.environ:
        print('exported PP_PATH=%s' % pp_path)
        os.environ['PP_PATH'] = pp_path

    path = os.environ['PATH']

    if pp_path not in path:
        os.environ['PATH'] = "%s:" + os.environ['PATH'] % pp_path
        print(os.environ['PATH'])
    else:
        print('PP_PATH is already in PATH')

    if found_argcomplete:
        subprocess.call('activate-global-python-argcomplete')
        subprocess.call('eval "$(register-python-argcomplete picpy.py)"', shell=True)
        subprocess.call('eval "$(register-python-argcomplete pp)"', shell=True)

    print('Consider adding the following lines in your bashrc:')
    print('export PP_PATH=%s' % pp_path)
    print('export PATH=$PP_PATH:$PATH')
    if found_argcomplete:
        print('eval "$(register-python-argcomplete picpy.py)"')
        print('eval "$(register-python-argcomplete pp)"')
        print('activate-global-python-argcomplete')

if __name__ == "__main__":
    main()






