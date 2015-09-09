
$useTarget = $ARGV[0];

if($useTarget==1){
	open(INTARGET, "<targetdistanceweights.dat");
}
open(ININITIAL,"<initialdistanceweights.dat");
open(INCURRENT,"<currentdistanceweights.dat");
open(INBESTCURRENT,"<currentbestdistanceweights.dat");

$counter=0;

%targetHash=();
%initialHash=();
%currentHash=();
%bestcurrentHash=();


#while(<INTARGET>){
while(<ININITIAL>){
	
	$init=$_;
	if($useTarget==1){
		$targ=<INTARGET>;
	}
	$curr=<INCURRENT>;
	$bestcurr=<INBESTCURRENT>;

	if($useTarget==1){
		$targ=~/([0-9.]+)\s([0-9.]+)/;
		$targkey=$1;
		$targval=$2;
	}
	$init=~/([0-9.]+)\s([0-9.]+)/;
	$initkey=$1;
	$initval=$2;
	$curr=~/([0-9.]+)\s([0-9.]+)/;
	$currkey=$1;
	$currval=$2;

	$bestcurr=~/([0-9.]+)\s([0-9.]+)/;
	$bestcurrkey=$1;
	$bestcurrval=$2;
	
	#if(check(\%targethash,$targkey)!=-1.0){
	if($useTarget==1){
		$targetHash{$targkey}=$targval;
	}
	#}
	#if(check(\%initalhash,$initkey)!=-1.0){
		$initialHash{$initkey}=$initval;
	#}
	#if(check(\%currenthash,$currkey)!=-1.0){
		$currentHash{$currkey}=$currval;
	#}
	$bestcurrentHash{$bestcurrkey}=$bestcurrval;

	#print $counter,"\t",$targval,"\t",$initval,"\t",$currval,"\n";
	#print $counter,"\t0.0\t",$initval,"\t",$currval,"\n";
	$counter=$counter+1;
}	
if($useTarget==1){
	close(INTARGET);
}
close(ININITIAL);
close(INCURRENT);
close(INBESTCURRENT);


if($useTarget==1){
	foreach $key (sort {$a<=>$b} keys %initialHash){
		print $key,"\t",$targetHash{$key},"\t",$initialHash{$key},"\t",$currentHash{$key},"\t",$bestcurrentHash{$key},"\n";
	}
}
else{
	foreach $key (sort {$a<=>$b} keys %initialHash){
		print $key,"\t0.0\t",$initialHash{$key},"\t",$currentHash{$key},"\t",$bestcurrentHash{$key},"\n";
	}
}


sub check {
    my ($href, $key) = @_;
    return (exists $href->{$key}) ? $href->{$key} : -1.0;
}

