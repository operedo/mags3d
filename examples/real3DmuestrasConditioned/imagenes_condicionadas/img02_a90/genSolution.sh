ex=$1
makedir=$2
#perl plotWeights.pl 0 > weights.dat
#perl plotVariograms.pl > variograms.dat
#awk < salida.out.log '{ print $6 }' > cost.dat
perl plotFiles.pl $ex 0 # gnuplot
perl plotPixelplt.pl $ex
if [ $makedir -eq 1 ]
then
	perl plotFiles.pl $ex 1 # octave
	#sleep 3
	mkdir $ex
	cp run.sh weights.dat variograms.dat cost.dat *image.m *image.dat *distanceweights.dat *variogram.dat salida.out.log salida.err.log $ex
	mv *_$ex.ps *_$ex.eps *_$ex.png  $ex
fi
