plot "density_r.dat" u 1:2 lc 4
set term png
set output "triangles_r.png"
set key off
set xlabel "r"
set ylabel "number of triangles"
set title "number of triangles vs distance to the center"
replot