#!/bin/bash
#To plot the cellDisplacement of a tip of the flap

gnuplot -persist > /dev/null 2>&1 << EOF
	set title "Y_Displacement of a flap vs. Time"
	set xlabel "Time / Iteration"
	set ylabel "Y_Displacement(m)"

	plot	"disp.txt" using 1:2 title 'Y Cell Displacement'
EOF

