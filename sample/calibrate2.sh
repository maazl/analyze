#!/bin/bash
echo Zero mode calibration consists of two parts.
echo In the first part the test impedance must be Z = 0 \(short cut\)
echo The second part requires Z = infinity \(open\).
echo You will get about 15 Sseconds when you are promted to change the setup.
echo It is essential that these steps are taken without reinitializing the sound device.
read -p "Press enter to start or Ctrl-C to abort."
./analyze \@calibrate2.cfg $@
echo The results have been saved to zero.dat.
echo You should copy this to a more meaningful name matching your setup.
gnuplot -persist setup.gp gpenv.gp zero.gp
