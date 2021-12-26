set grid
#set size ratio 1.4142 1,1.886
set y2range [.01:100]
set y2tics autofreq
set logscale xy
set logscale y2
plot [fmin:fmax] [1:10000] "data.dat" u 1:(-1000000000000./rref/2/pi/$9/$1) t "ESC [pF] <" w l \
  , 'data.dat' u 1:($8*rref/1000) t "ESR [kOhm] <" w l \
  , 'data.dat' u 1:(1000*$9*rref/$1/2/pi) t "ESL [mH] <" w l \
  , 'data.dat' u 1:(sqrt($9*$9+$8*$8)*rref/1000) t "|Z| [kOhm] <" w l \
  , 'data.dat' u 1:($8<0?-1/$8*$9*$1/1000:$8/$9*$1/1000) ax x1y2 t "fT [kHz] >" w l \
  , 'data.dat' u 1:($8/abs($9)) ax x1y2 t "tan delta [] >" w l
unset logscale
unset y2tics
