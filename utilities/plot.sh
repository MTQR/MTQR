#!/bin/bash

gnuplot -p <<- EOF

set terminal svg size 1800,1012
set output 'Estimate.svg'
set xlabel 'Lambda' font ",22"
set xlabel offset 0,-2.0
set xtics font ",15"
set xtics offset 0,-1,0
set xrange [-0.25:37]
set autoscale y
set ytics font ",18"
set logscale y
set grid
set key font ",28"
set datafile separator ","
plot "double.csv" u 1:2 with lines title "A-posteriori error (double precision)" lw 2 lt rgb "red", "quadruple.csv" u 1:2 with lines title "A-posteriori error (quadruple precision)" lw 2 lt rgb "blue", "epsilon.csv" u 1:2 with lines title "Machine espilon (double precision)" lw 1.5 lt rgb "black" 

EOF
	
	
