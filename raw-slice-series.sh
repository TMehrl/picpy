#!/bin/bash

if [ "$#" -eq 0 ]; then
  echo "No path specified, using ./DATA..." 
  DATA_PATH="DATA"
elif [ "$#" -eq 1 ]; then
  DATA_PATH=$1
else 
  echo "Illegal number of parameters!"    
fi

nohup raw-slice-series.py $DATA_PATH 1> rss.out 2> rss.err &