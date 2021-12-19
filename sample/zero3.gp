if (!exists("datafile")) datafile='zero.dat'
i={0,1}
unset y2tics
set format x "%.0s%c"
set grid
set multiplot
set size .5,.5
set origin 0,.5
plot [fmin:fmax] [.9:1.1] datafile u 1:10 t '|Cll|' w l, \
                 datafile u 1:16 t '|Crr|' w l
set origin 0,0
plot [fmin:fmax] [-10:10] datafile u 1:11 t 'arg Cll [°]' w l, \
                 datafile u 1:17 t 'arg Crr [°]' w l
set origin .5,.5
plot [fmin:fmax] [0:.3] datafile u 1:12 t '|Clr|' w l, \
                 datafile u 1:14 t '|Crl|' w l
set origin .5,0
#plot [fmin:fmax] [-180:180] datafile u 1:13 t 'arg Clr [°]' w p pt 0, \
#                 datafile u 1:15 t 'arg Crl [°]' w p pt 0
plot [fmin:fmax] [*:*] datafile u 1:18 t '|Uo|' w l, \
                 datafile u 1:24 t '|I0|' w l
unset multiplot

