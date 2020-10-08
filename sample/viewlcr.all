set grid
#set size noratio
set multiplot
set size .5,.5
#set y2range [.01:100]
#set y2tics autofreq
set logscale xy
#set logscale y2
set origin 0,.5
plot [fmin:fmax] [] "data.dat" u 1:6 t "|H|" w l
set origin 0,0
unset logscale
set logscale x
#plot [fmin:fmax] [] "data.dat" u 1:7 t "arg H" w l
plot [fmin:fmax] [] "data.dat" u 1:(($3-$5+360)-360*int(($3-$5+360)/360)) t "arg H" w l
set origin .5,.5
set logscale xy
plot [fmin:fmax] [] "data.dat" u 1:8 t "re H" w l
set origin .5,0
set logscale xy
plot [fmin:fmax] [] "data.dat" u 1:9 t "im H" w l
unset logscale
unset multiplot

