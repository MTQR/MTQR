#!/bin/bash

for sub_dir in ./*; do

	files=$(ls $sub_dir*.svg 2> /dev/null | wc -l)

	if [ "$files" != "1" ]; then

		gnuplot -p <<- EOF

		set terminal svg size 1800,1012
		set output '$sub_dir/Estimate$1.svg'
		set xlabel "\\(\\frac{x}{x_0}\\)" font ",22"
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
		plot "./estimate/EstimateN=$1.csv" u 1:2 with lines title "Exact estimate" lw 2 lt rgb "red", "./estimate/EnvelopedEstimateN=$1.csv" u 1:2 with lines title "Enveloped estimate" lw 2 lt rgb "blue", "./estimate/Epsilon.csv" u 1:2 with lines title "Machine d.p. espilon" lw 1.5 lt rgb "black"

		EOF

	fi
	
done