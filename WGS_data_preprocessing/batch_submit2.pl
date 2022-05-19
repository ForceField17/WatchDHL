#!/usr/bin/perl

#please run this job in a new clear fold
die "usage: perl maping.pl <Data_Index> <fastq-dir> <START> <END> \n" if @ARGV!= 4;

$index_table = $ARGV[0];
$root = $ARGV[1];


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


$start = $ARGV[2];
$end = $ARGV[3];

for($j=$start;$j<=$end;$j++){          
        `sh ./Run_BWA_mapping.sh $sample[$j][1] $root/$sample[$j][1]\.R1.fastq.gz $root/$sample[$j][1]\.R2.fastq.gz $sample[$j][1]`;
	`sh ./Run_mark_duplicates.sh $sample[$j][1]`;
	$completed_No = $j + 1;
        print "Job completed!!! No.$completed_No : $sample[$j][1]\t$sample[$j][0]\n\n";
}