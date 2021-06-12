#!/bin/bash
#To plot the pressure forces and Viscous forces

gnuplot -persist > /dev/null 2>&1 << EOF
	set title "Forces vs. Time"
	set xlabel "Time / Iteration"
	set ylabel "Force (N)"

	plot	"forces.txt" using 1:2 title 'Pressure Forces' with linespoints,\
			"forces.txt" using 1:3 title 'Viscous Forces' with linespoints
EOF