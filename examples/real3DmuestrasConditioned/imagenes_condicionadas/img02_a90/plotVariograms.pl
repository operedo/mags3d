
open(INTARGET,"<targetvariogram.dat");
#open(INTARGET,"<gamv_Cu_vertical.out");
#open(INTARGET,"<gamv_Cu_omnihoriz_nscore.out");
open(ININITIAL,"<initialvariogram.dat");
open(INCURRENT,"<currentvariogram.dat");
open(INBESTCURRENT,"<currentbestvariogram.dat");

$lines=<INTARGET>;
$lines=<ININITIAL>;
$lines=<INCURRENT>;
$lines=<INBESTCURRENT>;

$counter=0;

while($counter<$lines){
	$targ=<INTARGET>;
	$init=<ININITIAL>;
	$curr=<INCURRENT>;
	$bestcurr=<INBESTCURRENT>;

	
	$targ=~/([0-9.]+)\s+([0-9]+)\s+([0-9]+)/;
	$targval=$1;
	$targlag=$3;
	$init=~/([0-9.]+)\s+([0-9]+)\s+([0-9]+)/;
	$initval=$1;
        $initlag=$3;
	$curr=~/([0-9.]+)\s+([0-9]+)\s+([0-9]+)/;
	$currval=$1;
        $currlag=$3;
	$bestcurr=~/([0-9.]+)\s+([0-9]+)\s+([0-9]+)/;
	$bestcurrval=$1;
        $bestcurrlag=$3;
	
	print $targlag,"\t",$targval,"\t",$initlag,"\t",$initval,"\t",$currlag,"\t",$currval,"\t",$bestcurrlag,"\t",$bestcurrval,"\n";
	#print $counter,"\t",$targval,"\n";
	$counter=$counter+1;
}	
close(INTARGET);
close(ININITIAL);
close(INCURRENT);
close(INBESTCURRENT);
