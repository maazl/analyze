#!/bin/bash
# This file is called by analyze.
amixer -sq <mixer.in;

device="-D hw:0,0"
#device="-D sysdefault"
format="-r 48000 -f S16_LE -c 2 --disable-resample"
vumeter="-V stereo"
#recopt="--buffer-size 65536"

buffer="4M"

# examine current terminal emulation
term=`eval "$(xprop -notype -id "$WINDOWID" 32i '=$0' _NET_WM_PID)"; ps -o comm= -p "$_NET_WM_PID"`

alsaopt="$device $format $alsaopt"
recopt="$recopt $vumeter"
if [ -z "$term" ]; then
  if ! [ -p "play.pipe" ]; then
    mkfifo play.pipe;
    if [ -z "$buffer" ]; then
      bash -c "aplay $alsaopt -t raw play.pipe || true; sleep 1; rm play.pipe" &
    else
      bash -c "stdbuf -i$buffer aplay $alsaopt -t raw - <play.pipe || true; sleep 1; rm play.pipe" &
    fi;
  fi;
  if ! [ -p "rec.pipe" ]; then
    mkfifo rec.pipe;
    if [ -z "$buffer" ]; then
      bash -c "arecord $alsaopt $recopt -t raw >>rec.pipe || true; rm rec.pipe" &
    else
      bash -c "stdbuf -o$buffer arecord $alsaopt $recopt -t raw >>rec.pipe || true; rm rec.pipe" &
    fi;
  fi;
else
  if ! [ -p "play.pipe" ]; then
    mkfifo play.pipe;
    if [ -z "$buffer" ]; then
      $term -e "bash -c \"printf \\\"\033[8;5;80t\\\"; aplay $alsaopt -t raw play.pipe || true; sleep 1; rm play.pipe\"";
    else
      $term -e "bash -c \"printf \\\"\033[8;5;80t\\\"; stdbuf -i$buffer aplay $alsaopt -t raw - <play.pipe || true; sleep 1; rm play.pipe\"";
    fi;
  fi;
  if ! [ -p "rec.pipe" ]; then
    mkfifo rec.pipe;
    if [ -z "$buffer" ]; then
      $term -e "bash -c \"printf \\\"\033[8;5;80t\\\"; arecord $alsaopt $recopt -t raw >>rec.pipe || true; rm rec.pipe\"";
    else
      $term -e "bash -c \"printf \\\"\033[8;5;80t\\\"; stdbuf -o$buffer arecord $alsaopt $recopt -t raw >>rec.pipe || true; rm rec.pipe\"";
    fi;
  fi;
fi;
