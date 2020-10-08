#!/bin/bash
if [ ! "$1" ]; then
  echo "usage: $0 @<configfile>";
else
  nice -1 ./analyze $@ | nice -2 gnuplot gnuplot.in -;
fi;