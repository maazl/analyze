#!/bin/bash
# This file is called by analyze.
pulseaudio -k
amixer -c SB -sq <mixer.in;

device="-D hw:0,0"
#device="-D sysdefault"
#format="-r 192000 -f S16_LE -c 2 --disable-resample"
format="-r 48000 -f S16_LE -c 2 --disable-resample"
vumeter="-V stereo"
recopt="--buffer-size 1048576 -F 100000"

recpipe="/tmp/rec.pipe"
playpipe="/tmp/play.pipe"
#buffer="4M"

# examine current terminal emulation
term=`eval "$(xprop -notype -id "$WINDOWID" 32i '=$0' _NET_WM_PID)"; ps -o comm= -p "$_NET_WM_PID"`

alsaopt="$device $format $alsaopt"
recopt="$recopt $vumeter"
rm $recpipe $playpipe 2>/dev/null
if [ -z "$term" ]; then
  if ! [ -p "$playpipe" ]; then
    mkfifo $playpipe;
    if [ -z "$buffer" ]; then
      bash -c "aplay $alsaopt -t raw $playpipe || true; sleep 1; rm $playpipe" &
    else
      bash -c "stdbuf -i$buffer aplay $alsaopt -t raw - <$playpipe || true; sleep 1; rm $playpipe" &
    fi;
  fi;
  if ! [ -p "$recpipe" ]; then
    mkfifo $recpipe;
    if [ -z "$buffer" ]; then
      bash -c "arecord $alsaopt $recopt -t raw >>$recpipe || true; rm $recpipe" &
    else
      bash -c "stdbuf -o$buffer arecord $alsaopt $recopt -t raw >>$recpipe || true; rm $recpipe" &
    fi;
  fi;
else
  if ! [ -p "$playpipe" ]; then
    mkfifo $playpipe;
    if [ -z "$buffer" ]; then
      $term -e "bash -c \"printf \\\"\033[8;5;80t\\\"; aplay $alsaopt -t raw $playpipe || true; sleep 1; rm $playpipe\"";
    else
      $term -e "bash -c \"printf \\\"\033[8;5;80t\\\"; stdbuf -i$buffer aplay $alsaopt -t raw - <$playpipe || true; sleep 1; rm $playpipe\"";
    fi;
  fi;
  if ! [ -p "$recpipe" ]; then
    mkfifo $recpipe;
    if [ -z "$buffer" ]; then
      $term -e "bash -c \"printf \\\"\033[8;5;80t\\\"; arecord $alsaopt $recopt -t raw >>$recpipe || true; rm $recpipe\"";
    else
      $term -e "bash -c \"printf \\\"\033[8;5;80t\\\"; stdbuf -o$buffer arecord $alsaopt $recopt -t raw >>$recpipe || true; rm $recpipe\"";
    fi;
  fi;
fi;
