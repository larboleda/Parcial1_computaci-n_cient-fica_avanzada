plot "mass_profile.dat" u 2:(log($1)) lc 4
set term png
set output "cumulative_mass.png"
set key off
set xlabel "log(r)"
set ylabel "Cumulative Mass"
set title "Cumulative vs distance to the center"
replot