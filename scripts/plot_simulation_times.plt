#!/usr/bin/gnuplot --persist

set title "Test run provided in ALICIA-examples/ folder"
unset key
set grid

set xlabel "Shot time (s)"  
set ylabel "Clock time (minutes)" 

set xrange [46:]
plot "simulation_times.dat" u 1:($3/60E6) w l title "Simulation clock time"

