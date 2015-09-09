
open(IN,"<$ARGV[0]");
$mu=$ARGV[1];

$title=<IN>;
print $title;
$varnum=<IN>;
print $varnum;

for($var=0;$var<$varnum;$var++){
	$varname=<IN>;
	print $varname;
}

while(<IN>){
	$_=/([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)/;
	$x=$1;
	$y=$2;
	$z=$3;
	$data=$4;
	$data = $data + $mu;
	print $x," ",$y," ",$z," ",$data,"\n";
}

close(IN);
