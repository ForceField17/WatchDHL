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
$head = <FILE>;
chomp($head);
print "$head\tMin_CpG\tMin_CGC\tMin_WRCH\tMin_WRCY\tMin_WGCW\n";
$i=0;
while($line = <FILE>){
	chomp($line);
	$content[$i] = $line;
	@temp=split('\t',$line);
	$chr = $temp[1];
	$pos = $temp[2];

	$obj1 = $pos-500;
	$obj2 = $pos+500;	
	my $left = `samtools faidx /home/dsongad/software/my_lib/star_genome_d1_vd1_gtfv22/GRCh38.d1.vd1.fa $chr\:$obj1\-$obj2 | tail -n +2`;
	$left =~ s/\n//g;
	$left =~ s/[a-z]/[A-Z]/g;
	chomp($left);
	my @tmp=split('',$left);

	$min_CG = 500;
	$min_CGC = 500;
	$min_WRCH = 500;
	$min_WRCY = 500;
	$min_WGCW = 500;
	for($k=0;$k<=($#tmp-3);$k++){
		$motif_CG = $tmp[$k].$tmp[$k+1];
		$motif_CGC = $tmp[$k].$tmp[$k+1].$tmp[$k+2];
		$motif_WRCH = $tmp[$k].$tmp[$k+1].$tmp[$k+2].$tmp[$k+3];
	
                if( $motif_CG eq "CG" or $motif_CG eq "GC"){
                        if( $min_CG > abs($k-500+1) ){
				$min_CG = abs($k-500+1);
				if( $min_CG > abs($k-500) ){
                                	$min_CG = abs($k-500);
					if( $min_CG > abs($k-501) ){
						$min_CG = abs($k-501);
					}
				}
			}
                }

	        if( $motif_CGC eq "CGC" or $motif_CGC eq "GCG"){
                        if( $min_CGC > abs($k-500+2) ){
                                $min_CGC = abs($k-500+2);
                        }
                }

                if( $motif_WRCH =~ /[AT][AG]C[CTA]/ or $motif_WRCH =~ /[GAT]G[TC][AT]/){
                        if( $min_WRCH > abs($k-500+3) ){
                                $min_WRCH = abs($k-500+3);
                        }
                }

                if( $motif_WRCH =~ /[AT][AG]C[CT]/ or $motif_WRCH =~ /[GA]G[TC][AT]/){
                        if( $min_WRCY > abs($k-500+3) ){
                                $min_WRCY = abs($k-500+3);
                        }
                }

		if( $motif_WRCH =~ /[AT]GC[AT]/ or $motif_WRCH =~ /[AT]CG[AT]/){
                        if( $min_WGCW > abs($k-500+3) ){
                                $min_WGCW = abs($k-500+3);
                        }
                }

		
	}	
	
	 print "$line\t$min_CG\t$min_CGC\t$min_WRCH\t$min_WRCY\t$min_WGCW\n";	

	$i++;
}
close FILE;


