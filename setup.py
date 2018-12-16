#!/usr/bin/env python3

import os
import sys
import subprocess

def main():

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

    subprocess.call('activate-global-python-argcomplete')
    subprocess.call('eval "$(register-python-argcomplete picpy.py)"', shell=True)
    subprocess.call('eval "$(register-python-argcomplete pp)"', shell=True)


if __name__ == "__main__":
    main()






