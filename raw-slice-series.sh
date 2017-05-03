#!/bin/bash

nohup raw-slice-series.py "$@" 1> rss.out 2> rss.err &