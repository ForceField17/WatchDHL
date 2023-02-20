#!/usr/bin/perl

#please run this job in a new clear fold
die "usage: perl checking.pl <savi_table> \n" if @ARGV!= 1;

my %trans;
$trans{A} = 'T';
$trans{T} = 'A';
$trans{G} = 'C';
$trans{C} = 'G';
$trans{a} = 'T';
$trans{t} = 'A';
$trans{g} = 'C';
$trans{c} = 'G';
$trans{N} = 'N';
$trans{n} = 'N';

$target = $ARGV[0];
open FILE, $target;
#open FASTA,"> _$prefix\_fasta.fa";
#$head = <FILE>;
#chomp($head);
#print "$head\tforwardMutPattern\treverseMutPattern\tWRCY\n";
$i=0;
while($line = <FILE>){
	chomp($line);
	$content[$i] = $line;
	@temp=split('\t',$line);
	$chr = $temp[0];
	$pos = $temp[1]+1;

	$obj1 = $pos-2;
	$obj2 = $pos+2;

        $ref = $temp[3];
        $alt = $temp[4];

	if($temp[7] ne "SNV"){
		print "$line\tNA\n";
		next;
	}

        if( ($temp[3] eq "C" and $temp[4] eq "T") or ($temp[3] eq "G" and $temp[4] eq "A") ){
        
		my $left = `samtools faidx /home/dsongad/software/my_lib/star_genome_d1_vd1_gtfv22/GRCh38.d1.vd1.fa $chr\:$obj1\-$obj2 | tail -n +2`;
		$left =~ s/\n//g;
		$left =~ s/[a-z]/[A-Z]/g;
		chomp($left);
		$seqL[$i] = $left;
		my @tmp=split('',$left);
	
		$forward = $tmp[0].$tmp[1].$ref.$tmp[3];
		$reverse = $trans{$tmp[4]}.$trans{$tmp[3]}.$trans{$ref}.$trans{$tmp[1]};
	
		if( $forward =~ /^[AT][GA]C[TC]$/ or $reverse =~ /^[AT][GA]C[TC]$/){
			print "$line\tWRCY\n";
		}
		else{
			print "$line\tnonWRCY\n";
		}
	}
	else{
                print "$line\tNotC2T\n";
	}
}
close FILE;


