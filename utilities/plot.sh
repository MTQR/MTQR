#!/bin/bash

for sub_dir in ../test/*/; do

	files=$(ls $sub_dir*.svg 2> /dev/null | wc -l)

	if [ "$files" != "1" ]; then

		gnuplot -p <<- EOF

		set terminal svg size 1800,1012
		set output '$sub_dir/error_estimate_plot(n=$1).svg'
		set xlabel "lambda"
		set xrange [$2:$3]
		set autoscale y
		set ylabel "Error estimate"
		set logscale y
		set grid
		set datafile separator ","
		plot "../test/output/ExactError.csv" u 1:2 with lines title "Exact estimate" lw 1.5 lt rgb "red", "../test/output/ExactErrorEnv.csv" u 1:2 with lines title "Enveloped estimate" lw 1.5 lt rgb "blue", "../test/output/Epsilon.csv" u 1:2 with lines title "Precision (double-format)" lw 1.5 lt rgb "black" 

		EOF

	fi
	
done
	
	
