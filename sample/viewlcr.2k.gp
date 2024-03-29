set grid
#set size ratio 1.4142 1,1.886
set y2range [.01:100]
set y2tics autofreq
set logscale xy
set logscale y2
plot [fmin:fmax] [.1:1000] "data.dat" u 1:(-1000000000/rref/2/pi/$9/$1) t "C [nF] <" w l \
  , 'data.dat' u 1:($8*rref/1000) ax x1y2 t "R [kOhm] >" w l \
  , 'data.dat' u 1:($9*rref/$1/2/pi) ax x1y2 t "L [H] >" w l \
  , 'data.dat' u 1:($8<0?-1/$8*$9*$1/1000:$8/$9*$1/1000) ax x1y2 t "fT [kHz] >" w l \
  , 'data.dat' u 1:($8/abs($9)) ax x1y2 t "tan delta [] >" w l \
  , 'data.dat' u 1:(sqrt($9*$9+$8*$8)*rref/1000) t "|Z| [kOhm] <" w l
unset logscale
unset y2tics
