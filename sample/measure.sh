#!/bin/bash
if [ ! "$1" ]; then
  echo "usage: $0 @<configfile>";
else
  nice -5 ../analyze $@ | nice -10 gnuplot setup.gp -;
fi;