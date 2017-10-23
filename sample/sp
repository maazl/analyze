set logscale x
set grid
set size 1,1
set origin 0,0
set multiplot
set size .5,.5
set origin 0,.5
dist=1.02/320.
plot [20:20000] 'log' u 1:(20.*log($6)/log(10)) t 'abs' w l
set origin 0,0
plot [20:20000] 'log' u 1:($7-360.*dist*$1-360*floor($7/360.+.5-dist*$1)) t 'arg' w p pt 0
set origin .5,.5
plot [20:20000] 'log' u 1:(20*log($2/$12)/log(10)) t 'SNRref' w l
set origin .5,0
plot [20:20000] 'log' u 1:(20*log($4/$13)/log(10)) t 'SNRmic' w l
set nomultiplot
