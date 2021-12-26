set grid
#set size ratio 1.4142 1,1.886
set y2range [.01:100]
set y2tics autofreq
set logscale xy
set logscale y2
plot [fmin:fmax] [.1:1000] "data.dat" u 1:(-1000/rref/2/pi/$9/$1) ax x1y2 t "ESC [mF] >" w l \
  , 'data.dat' u 1:(1000*$8*rref) t "ESR [mOhm] <" w l \
  , 'data.dat' u 1:($9*rref/$1/2/pi*1000000) t "ESL [µH] >" ax x1y2 w l \
  , 'data.dat' u 1:(1000*$6*rref) t "|Z| [mOhm] <" w l \
  , 'data.dat' u 1:($9<0?-1*$9/$8*$1/1000:$8/$9*$1/1000) t "fT [kHz] <" w l \
  , 'data.dat' u 1:($8/abs($9)) ax x1y2 t "|tan delta| [] >" w l
unset logscale
unset y2tics
