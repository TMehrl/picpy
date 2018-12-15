#!/usr/bin/env python3

import sys
import argparse
import argcomplete


def main(**args):
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('positional', choices=['g3d', 'raw2hist', 'rss', 'rss-plot'])
    parser.add_argument('--optional', choices=['foo1', 'foo2', 'bar'])
    argcomplete.autocomplete(parser)
    args = parser.parse_args()
    main(**vars(args))
