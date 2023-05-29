#!/bin/bash

for sub_dir in ./*; do

	files=$(ls $sub_dir*.svg 2> /dev/null | wc -l)

	if [ "$files" != "1" ]; then

		gnuplot -p <<- EOF

		set terminal svg noenhanced size 1800,1012
		set output '$sub_dir/Fig2.svg'
		set xlabel "λ" font ",22"
		set xlabel offset 0,-2.0
		set xtics font ",12"
		set xtics offset 0,-1,0
		set xrange [$2:$3]
		set autoscale y
		set ytics ("$10^{-25}$" 1E-25, "$10^{-20}$" 1E-20, "$10^{-15}$" 1E-15, "$10^{-10}$" 1E-10, "$10^{-5}$" 1E-5, "$1$" 1) 
		set ytics font ",12"
		set logscale y
		set grid
		set key font ",20"
		set key right top
		set datafile separator ","
		plot "./estimate/EstimateN=$1.csv" u 1:2 with lines title "Exact estimate" lw 2 lt rgb "red", "./estimate/EnvelopedEstimateN=$1.csv" u 1:2 with lines title "Enveloped estimate" lw 2 lt rgb "blue", "./estimate/Epsilon.csv" u 1:2 with lines title "Machine d.p. espilon" lw 1.5 lt rgb "black"

		EOF

	fi
	
done