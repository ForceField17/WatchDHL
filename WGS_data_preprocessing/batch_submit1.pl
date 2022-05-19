#!/usr/bin/perl

#please run this job in a new clear fold
die "usage: perl checking.pl <Data_Index> <start> <end> \n" if @ARGV!= 3;

$index_table = $ARGV[0];



open INDEX, $index_table;
$i = 0;
$id = 0;
while($line = <INDEX>){
	chomp($line);
	my @temp=split('\t',$line);
	$sample[$i][0] = $temp[1];
	$sample[$i][1] = $temp[2];
	$i++;
}
close INDEX;


$start = $ARGV[1];
$end = $ARGV[2];

for($j=$start;$j<=$end;$j++){          #run the first 100 samples as a try
        `sh ./Run_fastp_readsFilter.sh $sample[$j][0] $sample[$j][1]`;
	$completed_No = $j;
        print "Job completed!!! No.$completed_No : $sample[$j][1]\t$sample[$j][0]\n\n";
}


