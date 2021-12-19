#!/bin/bash
echo VERIFY ZERO MODE CALIBRATION
echo Reads the last zero.dat and shows the residual difference.
echo In the first part the test impedance must be Z = 1
echo In the second part the test impedance must be Z = infinity \(open\)
echo And the last part requires Z = 0 \(short cut\).
echo You will get about 15 Sseconds when you are promted to change the setup.
echo It is essential that these steps are taken without reinitializing the sound device.
read -p "Press enter to start or Ctrl-C to abort."
./analyze \@calibrate3.cfg zr=zero.dat zg=zeroD.dat "$@"
gnuplot -persist setup.gp gpenv.gp zeroD.gp
