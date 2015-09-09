ex=$1
makedir=$2
perl plotWeights.pl 0 > weights.dat
perl plotVariograms.pl > variograms.dat
awk < salida.out.log '{ print $6 }' > cost.dat
perl plotFiles.pl $ex 0 # gnuplot
if [ $makedir -eq 1 ]
then
	#perl plotFiles.pl $ex 1 # octave
	#sleep 3
	mkdir $ex
	cp weights.dat variograms.dat cost.dat *image.m *distanceweights.dat *variogram.dat salida.out.log salida.err.log $ex
	mv *_$ex.png  $ex
fi
