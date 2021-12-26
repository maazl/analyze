set grid
#set size ratio 1.4142 1,1.886
#set y2range [.01:100]
set logscale x
unset logscale y
set autoscale xy
set y2tics autofreq
set autoscale y2
#set logscale y2
plot "<awk '$12==0 { print $0 }' data.dat" u 1:(20*log10($6)) w l t "|Hl| [dB] <", \
     "<awk '$12==1 { print $0 }' data.dat" u 1:(20*log10($6)) w l t "|Hr| [dB] <", \
     "<awk '$12==0 { print $0 }' data.dat" u 1:(1000*$11) w l ax x1y2 t "|delay l| [ms] >", \
     "<awk '$12==1 { print $0 }' data.dat" u 1:(1000*$11) w l ax x1y2 t "|delay r| [ms] >"
#unset logscale
#unset y2tics
