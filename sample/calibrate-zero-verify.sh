#!/bin/bash
echo VERIFY ZERO MODE CALIBRATION
echo Reads the last zero.dat and shows the residual difference.
echo In the first part the test impedance must be Z = 0 \(short cut\).
echo The second part requires Z = infinity \(open\).
echo You will get about 15 Sseconds when you are promted to change the setup.
echo It is essential that these steps are taken without reinitializing the sound device.
read -p "Press enter to start or Ctrl-C to abort."
./analyze \@calibrate.cfg zr=zero.dat zg=zeroD.dat "$@"
gnuplot -persist gpenv zeroD
