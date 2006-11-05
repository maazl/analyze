/**/
file = ARG(1)
IF file = '' THEN file = 'dataE.dat'
pd = LASTPOS('.',file)
pc = LASTPOS(':',file)
pb = LASTPOS('\',file)
IF pd <= pc | pd <= pb THEN pd = LENGTH(file)
ofile = SUBSTR(file, 1, pd-1)'2'SUBSTR(file, pd)
'perl intgrid.pl 10 .022091685 <"'file'" >"'ofile'"'
/* perl intgrid.pl 10 .029731547 <dataE.dat >dataE2.dat
 perl intgrid.pl 10 .029731547 <dataE.bak >dataE2.bak*/