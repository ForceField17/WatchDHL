#!/usr/bin/perl

#please run this job in a new clear fold
die "usage: perl checking.pl <savi_table> \n" if @ARGV!= 1;


$target = $ARGV[0];
open FILE, $target;
#open FASTA,"> _$prefix\_fasta.fa";
#$head = <FILE>;
#chomp($head);
#print "$head\n";
$i=0;
while($line = <FILE>){
	chomp($line);
	#$content[$i] = $line;
	@temp=split('\t',$line);
	#$gene = $temp[16];
	#$AA = $temp[17];
	#$Eff = $temp[22];
	$Imp = $temp[17];
		
	my @tmp = split(',',$Imp);
#	if($tmp[1] ne 'MODERATE' and $tmp[1] ne 'HIGH'){
		$line =~ s/\+[^\t]+\t/\t/g;
#	}
	print "$line\n";
}
close FILE;


