set term png;
plot 'weights.dat' u 1:2 t 'Target' w lp, 'weights.dat' u 1:3 t 'Initial' w lp, 'weights.dat' u 1:4 t 'Current' w lp, 'weights.dat' u 1:5 t 'Best' w lp;
set output "weights_sol01.png";
replot;
plot 'variograms.dat' u 1:2 t 'Target' w lp, 'variograms.dat' u 1:3 t 'Initial' w lp, 'variograms.dat' u 1:4 t 'Current' w lp, 'variograms.dat' u 1:5 t 'Best' w lp;
set output "variograms_sol01.png";
replot;
set logscale y;
plot 'cost.dat' u 1;
set output "cost_function_sol01.png";
replot;

