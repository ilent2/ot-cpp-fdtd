#!/usr/bin/gnuplot

set term png size 1024, 1024
F(x, y, z) = abs(y)

unset tics
unset colorbox

set size ratio -1

set pm3d map

set output "output_0.png"
splot "plane_0.dat" u 2:3:(F($4, $5, $6))

set output "output_1.png"
splot "plane_1.dat" u 2:3:(F($4, $5, $6))

