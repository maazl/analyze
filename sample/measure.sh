#!/bin/bash
if [ ! "$1" ]; then
  echo "usage: $0 @<configfile>";
else
  nice -1 ./analyze $@ | gnuplot gnuplot.in -;
fi;