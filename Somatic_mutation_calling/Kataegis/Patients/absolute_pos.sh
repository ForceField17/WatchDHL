#!/usr/bin/perl

#please run this job in a new clear fold

$index_table = $ARGV[0];


open XXX,"/home/dsongad/IGCT/run_savi/all_germSNP2/Annotation_for_allSNP/step6_rmSecondaryAnno/MappingQ0/lib/coordinate.txt";
<XXX>;
$i=0;
while($line = <XXX>){
        chomp($line);
        my @yyy=split('\t',$line);
        $sample[$i][0] = $yyy[1];
        $sample[$i][1] = $yyy[3];
#	print "$sample[$i][1]\n";
	$i++;
}
close XXX;







open INDEX, $index_table;

$i = 0;
print "Chr\tstart\tCase\tref\talt\ttype\tk27\tme3\tboth2\tse\tAID\tlabel\thistone\tdis\tpos\n";
#print "$upper_line\n";

$count = 0;
$line = <INDEX>;
        chomp($line);
        my @temp=split('\t',$line);
        $chr = $temp[0];
        $start = $temp[1];
print "$line\t$start\n";
$mark = 0;
$i=0;
while($line = <INDEX>){
        chomp($line);
        my @yyy=split('\t',$line);
	$chr_new = $yyy[0];
	$start_new = $yyy[1];
	if($chr_new eq $chr){
		$pos = $mark + $start_new;
		print "$line\t$pos\n";
	}
	else{	
		for($i=0;$i<=$#sample;$i++){
			#print "$sample[$i][0]\t$sample[$i][1]\n";
			if($sample[$i][0] eq $chr){
				$mark = $sample[$i][1];
				last;
			}
		}
		#print "$mark\n";
		$pos = $mark + $start_new;
		print "$line\t$pos\n";
	}
	$start = $start_new;
	$chr = $chr_new;
}
close INDEX;

