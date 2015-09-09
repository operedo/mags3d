
$ex=$ARGV[0];

$mode=$ARGV[1];

if($mode==0){
open(GNUPLOT,">plotter.gnuplot");
$text=<<EOF;
set term png;
plot 'weights.dat' u 1:2 t 'Target' w lp, 'weights.dat' u 1:3 t 'Initial' w lp, 'weights.dat' u 1:4 t 'Current' w lp, 'weights.dat' u 1:5 t 'Best' w lp;
set output "weights_$ex.png";
replot;
plot 'variograms.dat' u 1:2 t 'Target' w lp, 'variograms.dat' u 1:3 t 'Initial' w lp, 'variograms.dat' u 1:4 t 'Current' w lp, 'variograms.dat' u 1:5 t 'Best' w lp;
set output "variograms_$ex.png";
replot;
set logscale y;
plot 'cost.dat' u 1;
set output "cost_function_$ex.png";
replot;
EOF
print GNUPLOT $text,"\n";
close(GNUPLOT);


$ret=`echo "load 'plotter.gnuplot'" | gnuplot`;
}

#$imageDim="reshape(image(:,:,1),100,100)";
#$imageDim="reshape(image(:,:,1),50,50)";
$imageDim="image(:,:,1)";

if($mode==1){
open(OCTAVE,">plotter_octave.m");

$text=<<EOF;
targetimage;
imagesc($imageDim);
print -dpng targetimage_$ex.png;
initialimage;
imagesc($imageDim);
print -dpng initialimage_$ex.png;
currentimage;
imagesc($imageDim);
print -dpng currentimage_$ex.png;
EOF

#$text=<<EOF;
#initialimage;
#imagesc($imageDim);
#print -dpng initialimage_$ex.png;
#currentimage;
#imagesc($imageDim);
#print -dpng currentimage_$ex.png;
#EOF

print OCTAVE $text,"\n";
close(OCTAVE);

$ret=`octave --eval 'plotter_octave'`;
}

#$ret=`mkdir $ex`;
#$ret=`mv *_$ex.png weights.dat variograms.dat cost.dat $ex`;
#$ret=`mv *image.m *distanceweights.dat *variogram.dat salida.out.log salida.err.log $ex`;



