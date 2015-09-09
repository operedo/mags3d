# set term png;
#set terminal postscript eps enhanced color font 'Helvetica,36'
set terminal postscript eps enhanced color font 'Helvetica,26'
set key right top;
set xlabel 'Radial distance';
set ylabel 'Weight value';
#set xtics 0,100,300;
#plot 'weights.dat' u 1:3 t 'Initial' w lp lc rgb "green" lw 5 ps 2, 'weights.dat' u 1:5 t 'Best' w p lc rgb "blue" lw 5 ps 0.5;
plot 'weights.dat' u 1:3 t 'Initial' w lp lc rgb "green", 'weights.dat' u 1:5 t 'Best' w lp lc rgb "blue";
set output "weights_img01_a150.eps";
replot;
unset output;
set key right bottom;
set xlabel 'Lag distance';
set ylabel 'Semivariogram value';
#set xtics 0,35,140;
#plot 'variograms.dat' u 1:2 t 'Target' w lp lw 5 ps 2, 'variograms.dat' u 3:4 t 'Initial' w lp lw 5 ps 2, 'variograms.dat' u 5:6 t 'Best' w lp lw 5 ps 2;
plot 'variograms.dat' u 1:2 t 'Target' w lp, 'variograms.dat' u 1:3 t 'Initial' w lp, 'variograms.dat' u 1:4 t 'Unconditioned' w lp, 'variograms.dat' u 1:7 t 'Conditioned' w lp;
set output "variograms_img01_a150.eps";
replot;
unset output;
set logscale y;
set key right top;
set xlabel 'Iterations';
set ylabel 'Current cost / Initial cost';
#set xtics 0,6000,12000;
#plot 'cost.dat' u 1 t 'Cost function' lw 5 ps 2;
plot 'cost.dat' u 1 t 'Cost function';
set output "cost_function_img01_a150.eps";
replot;
unset output;

