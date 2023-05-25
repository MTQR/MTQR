#!/bin/bash

for sub_dir in ./*; do

	files=$(ls $sub_dir*.svg 2> /dev/null | wc -l)

	if [ "$files" != "1" ]; then

		gnuplot -p <<- EOF

		set terminal svg size 1800,1012
		set output '$sub_dir/PlotN=$1.svg'
		set xlabel 'λ' font ",22"
		set xlabel offset 0,-2.0
		set xtics font ",15"
		set xtics offset 0,-1,0
		set xrange [$2:$3]
		set autoscale y
		set ytics font ",18"
		set logscale y
		set grid
		set key font ",28"
		set key right bottom
		set datafile separator ","
		plot "./estimate/EstimateN=$1.csv" u 1:2 with lines title "R_{N=$1}" lw 2 lt rgb "red", plot "./estimate/EstimateN=$1.csv" u 1:2 with lines title "R_{N=$1}" lw 2 lt rgb "red", "./estimate/Epsilon.csv" u 1:2 with lines title "Machine d.p. espilon" lw 1.5 lt rgb "black"

		EOF

	fi
	
done