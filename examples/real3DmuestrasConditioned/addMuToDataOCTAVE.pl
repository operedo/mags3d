
open(IN,"<$ARGV[0]");
$mu=$ARGV[1];

while(<IN>){
	$_=/(.+)=\s+([-0-9.]+);/;
	$image=$1;
	$data=$2;
	$data = $data + $mu;
	print $image,"= ",$data,";\n";
}

close(IN);
