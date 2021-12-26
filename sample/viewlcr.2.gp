set grid
#set size ratio 1 1,1
set y2range [.1:1000]
set y2tics autofreq
set logscale xy
set logscale y2
plot [fmin:fmax] [.001:10] "data.dat" u 1:(-1000/rref/2/pi/$9/$1) t "ESC [mF] <" w l \
  , 'data.dat' u 1:($8*rref) t "ESR [Ohm] <" w l \
  , 'data.dat' u 1:($9*rref/$1/2/pi*1000) t "ESL [mH] <" w l \
  , 'data.dat' u 1:($9<0?-1/$8*$9*$1/1000:$8/$9*$1/1000) ax x1y2 t "fT [kHz] >" w l \
  , 'data.dat' u 1:(abs($9)/$8) ax x1y2 t "Q [] >" w l \
  , 'data.dat' u 1:($6*rref) t "|Z| [Ohm] <" w l
#unset logscale
#unset y2tics
