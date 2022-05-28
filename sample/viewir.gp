unset logscale x
unset logscale y
plot 'ir.dat' u 1:2 w l t "left",\
	'ir.dat' u 1:3 w l t "right"