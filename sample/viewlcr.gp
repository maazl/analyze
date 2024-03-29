set grid
#set size ratio 1.4142 1,1.886
set y2range [.01:100]
set y2tics autofreq
set logscale xy
set logscale y2
plot [fmin:fmax] [.1:1000] "data.dat" u 1:(-1000000/rref/2/pi/$9/$1) t "ESC [µF] <" w l \
  , 'data.dat' u 1:($8*rref) t "ESR [Ohm] <" w l \
  , 'data.dat' u 1:($9*rref/$1/2/pi*1000) t "ESL [mH] >" ax x1y2 w l \
  , 'data.dat' u 1:($6*rref) t "|Z| [Ohm] <" w l \
  , 'data.dat' u 1:($9<0 ? -1/$8*$9*$1/1000 : $8/$9*$1/1000) t "fT [kHz] <" w l \
  , 'data.dat' u 1:($8/abs($9)) t "tan delta [] >" ax x1y2 w l
#unset logscale
#unset y2tics
