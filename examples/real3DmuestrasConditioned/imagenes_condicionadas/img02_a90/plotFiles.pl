
$ex=$ARGV[0];

$mode=$ARGV[1];

if($mode==0){
open(GNUPLOT,">plotter.gnuplot");
$text=<<EOF;
# set term png;
#set terminal postscript eps enhanced color font 'Helvetica,36'
set terminal postscript eps enhanced color font 'Helvetica,26'
set key right top;
set xlabel 'Radial distance';
set ylabel 'Weight value';
#set xtics 0,100,300;
#plot 'weights.dat' u 1:3 t 'Initial' w lp lc rgb "green" lw 5 ps 2, 'weights.dat' u 1:5 t 'Best' w p lc rgb "blue" lw 5 ps 0.5;
plot 'weights.dat' u 1:3 t 'Initial' w lp lc rgb "green", 'weights.dat' u 1:5 t 'Best' w lp lc rgb "blue";
set output "weights_$ex.eps";
replot;
unset output;
set key right bottom;
set xlabel 'Lag distance';
set ylabel 'Semivariogram value';
#set xtics 0,35,140;
#plot 'variograms.dat' u 1:2 t 'Target' w lp lw 5 ps 2, 'variograms.dat' u 3:4 t 'Initial' w lp lw 5 ps 2, 'variograms.dat' u 5:6 t 'Best' w lp lw 5 ps 2;
plot 'variograms.dat' u 1:2 t 'Target' w lp, 'variograms.dat' u 1:3 t 'Initial' w lp, 'variograms.dat' u 1:4 t 'Unconditioned' w lp, 'variograms.dat' u 1:7 t 'Conditioned' w lp;
set output "variograms_$ex.eps";
replot;
unset output;
set logscale y;
set key right top;
set xlabel 'Iterations';
set ylabel 'Current cost / Initial cost';
#set xtics 0,6000,12000;
#plot 'cost.dat' u 1 t 'Cost function' lw 5 ps 2;
plot 'cost.dat' u 1 t 'Cost function';
set output "cost_function_$ex.eps";
replot;
unset output;
EOF
print GNUPLOT $text,"\n";
close(GNUPLOT);


$ret=`echo "load 'plotter.gnuplot'" | gnuplot`;
}

#$xtick="[0, 20,40, 60, 80, 100]";
#$xticklabel="[\"0\";\"200\";\"400\";\"600\";\"800\";\"1000\"]";
#$ytick="[0, 20,40, 60, 80, 100]";
#$yticklabel="[\"0\";\"200\";\"400\";\"600\";\"800\";\"1000\"]";

if($mode==1){
open(OCTAVE,">plotter_octave.m");

$cmin=-2.0;
$cmax=4.0;

for($slice=1;$slice<=12;$slice++){

$im="reshape(image(:,:,$slice),40,60)";

$text=<<EOF;
currentimage_scaled;
figure;
imagesc($im);
#ax=gca();
#set (ax, "xtick", $xtick);
#set (ax, "xticklabel", $xticklabel); 
#set (ax, "ytick", $ytick);
#set (ax, "yticklabel", $yticklabel); 
title("Unconditional image slice=$slice (scaled)","fontsize",26);
xlabel("X","fontsize",20);
ylabel("Y","fontsize",20);
caxis([ $cmin , $cmax ]);
colorbar();
axis("xy");
print -dpng currentimage_scaled_$slice-slice_$ex.png;
conditionedimage_scaled;
figure;
imagesc($im);
#ax=gca();
#set (ax, "xtick", $xtick);
#set (ax, "xticklabel", $xticklabel); 
#set (ax, "ytick", $ytick);
#set (ax, "yticklabel", $yticklabel); 
title("Conditional image slice=$slice (scaled)","fontsize",26);
xlabel("X","fontsize",20);
ylabel("Y","fontsize",20);
caxis([ $cmin , $cmax ]);
colorbar();
axis("xy");
print -dpng conditionedimage_scaled_$slice-slice_$ex.png;
EOF
print OCTAVE $text,"\n";

}

#$slice=6;
#
#$im="reshape(image(:,:,$slice),40,60)";
#
#
#$text=<<EOF;
#currentimage_scaled;
#figure;
#imagesc($im);
##ax=gca();
##set (ax, "xtick", $xtick);
##set (ax, "xticklabel", $xticklabel); 
##set (ax, "ytick", $ytick);
##set (ax, "yticklabel", $yticklabel); 
#title("Unconditional image slice=$slice (scaled)","fontsize",26);
#xlabel("X","fontsize",20);
#ylabel("Y","fontsize",20);
#caxis([ $cmin , $cmax ]);
#colorbar();
#axis("xy");
#print -dpng currentimage_scaled_$slice-slice_$ex.png;
#conditionedimage_scaled;
#figure;
#imagesc($im);
##ax=gca();
##set (ax, "xtick", $xtick);
##set (ax, "xticklabel", $xticklabel); 
##set (ax, "ytick", $ytick);
##set (ax, "yticklabel", $yticklabel); 
#title("Conditional image slice=$slice (scaled)","fontsize",26);
#xlabel("X","fontsize",20);
#ylabel("Y","fontsize",20);
#caxis([ $cmin , $cmax ]);
#colorbar();
#axis("xy");
#print -dpng conditionedimage_scaled_$slice-slice_$ex.png;
#EOF
#print OCTAVE $text,"\n";
#
#
#
#$slice=12;
#
#$im="reshape(image(:,:,$slice),40,60)";
#
#$text=<<EOF;
#currentimage_scaled;
#figure;
#imagesc($im);
##ax=gca();
##set (ax, "xtick", $xtick);
##set (ax, "xticklabel", $xticklabel); 
##set (ax, "ytick", $ytick);
##set (ax, "yticklabel", $yticklabel); 
#title("Unconditional image slice=$slice (scaled)","fontsize",26);
#xlabel("X","fontsize",20);
#ylabel("Y","fontsize",20);
#caxis([ $cmin , $cmax ]);
#colorbar();
#axis("xy");
#print -dpng currentimage_scaled_$slice-slice_$ex.png;
#conditionedimage_scaled;
#figure;
#imagesc($im);
##ax=gca();
##set (ax, "xtick", $xtick);
##set (ax, "xticklabel", $xticklabel); 
##set (ax, "ytick", $ytick);
##set (ax, "yticklabel", $yticklabel); 
#title("Conditional image slice=$slice (scaled)","fontsize",26);
#xlabel("X","fontsize",20);
#ylabel("Y","fontsize",20);
#caxis([ $cmin , $cmax ]);
#colorbar();
#axis("xy");
#print -dpng conditionedimage_scaled_$slice-slice_$ex.png;
#EOF
#print OCTAVE $text,"\n";


close(OCTAVE);

$ret=`octave --eval 'plotter_octave'`;
}

#$ret=`mkdir $ex`;
#$ret=`mv *_$ex.png weights.dat variograms.dat cost.dat $ex`;
#$ret=`mv *image.m *distanceweights.dat *variogram.dat salida.out.log salida.err.log $ex`;



