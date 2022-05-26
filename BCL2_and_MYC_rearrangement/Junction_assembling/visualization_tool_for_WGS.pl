#!/usr/bin/perl
#programmed by SONG,DONG 2019.12, and updated in 2022.05 
#please run this job in a new clean fold
die "usage: perl fusion_visual.pl <report_file_1_line> <reference_genome_dir> <ReadsPool>\n" if @ARGV!= 3;

use List::Util qw(min max);
#use warnings;

$gjf = $ARGV[0];
$genome = $ARGV[1];
$pool = $ARGV[2];

@sampleName = split('\.',$pool);
print "\n";
open INDEX, $gjf;

use Term::ANSIColor qw(:constants);

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
$trans{_} = '-';



@negChain5 = ();
@posChain5 = ();


my $line = <INDEX>;
chomp($line);
my @temp=split('\s+',$line);
$geneLeft = $temp[2];
$geneRight = $temp[3];

@xx1 = split(':',$temp[6]);
@xx2 = split(':',$temp[7]); 

$chromosome5 = "chr".$xx1[0];
$chromosome3 = "chr".$xx2[0];
$breakpoint5 = $xx1[1];
$breakpoint3 = $xx2[1];

@yy1 = split('/',$temp[4]);
@yy2 = split('/',$temp[5]);

$geneDirL = $yy1[0];
$geneDirR = $yy2[0];

$chain5 = $yy1[1];
$chain3 = $yy2[1];
$reads_file = $pool;

open R1, "$reads_file";
$i=0;
while($line=<R1>){
	chomp($line);
	$data1[$i][0] = $line;
	$data1[$i][1] = -500;
	$i++;
}
close R1;
$i=0;
for($j=0;$j<=$#data1;$j++){
	@vvv = split("",$data1[$j][0]);
	$rrr="";
	for($k=$#vvv;$k>=0;$k--){
		$rrr .= $trans{$vvv[$k]};
	}
	$data2[$i][0] = $rrr;
	$data2[$i][1] = -500;
	$i++;
}

$i=0;
for($n=0;$n<=$#data1;$n++){
	$data[$i][0] = $data1[$n][0];
	$data[$i][1] = $data1[$n][1];
	$i++;
	$data[$i][0] = $data2[$n][0];
        $data[$i][1] = $data2[$n][1];
	$i++;
}


$insertion = $temp[9];
close INDEX;
print "$temp[0]\t$chromosome5:$breakpoint5 $geneLeft（$geneDirL) --- $geneRight ($geneDirR) $chromosome3:$breakpoint3\n";


$page1 = int($breakpoint5/50) - 1;
$page2 = $page1 + 6;   #readin 50 * 6 bp reference sequence
#print "$page1\t$page2\n";
`sed -n "$page1,$page2\p" $genome/$chromosome5\.fa > ./tmp.5.$ARGV[0]`;
open GENOME, "./tmp.5.$ARGV[0]";


$line_before1 = <GENOME>;
$line_before2 = <GENOME>;
$line_before3 = <GENOME>;
$line=<GENOME>;
$line_after1=<GENOME>;
$line_after2=<GENOME>;
$line_after3=<GENOME>;

chomp($line_before1);
chomp($line_before2);
chomp($line_before3);
chomp($line);
chomp($line_after1);
chomp($line_after2);
chomp($line_after3);

my $line_local = uc($line_before1.$line_before2.$line_before3.$line.$line_after1.$line_after2.$line_after3);
#print "$line_local\n";
my @temp=split('',$line_local);
$position = $breakpoint5 % 50 -1 ;
$coord = $position + length($line_before1) + length($line_before2) + length($line_before3);
$a = 0;
for($z=-140;$z<=140;$z++){
        $seq5[$a] = $temp[$coord+$z];
	$a++;
} 
close GENOME;
#`rm ./tmp.5.$ARGV[0]`;

$chime = "";
if($chain5 eq '+'){
	printf " Left:%15s + 5' --- ",$geneLeft;
	for($i=(2+$insertion);$i<=140;$i++){
		print BOLD MAGENTA "$seq5[$i]", RESET;
		$chime .= $seq5[$i]; 
	}
        for($i=141;$i<=$#seq5;$i++){
                print  "$seq5[$i]";
        }
	print " --- 3'\n";
}
elsif($chain5 eq '-'){
	printf " Left:%15s - 5' --- ",$geneLeft;
        for($i=($#seq5-1-$insertion);$i>140;$i--){
                print BOLD MAGENTA "$trans{$seq5[$i]}", RESET;
		$chime .= $trans{$seq5[$i]};
        }
	for($i=140;$i>=1;$i--){
		print "$trans{$seq5[$i]}";
	}
        print " --- 3'\n";
}
else{
	print "Warning!!! Wrong chain signal\n";
}




$page1 = int($breakpoint3/50) - 1;
$page2 = $page1 + 6;
`sed -n "$page1,$page2\p" $genome/$chromosome3\.fa > ./tmp.3.$ARGV[0]`;
open GENOME, "./tmp.3.$ARGV[0]";


$line_before1 = <GENOME>;
$line_before2 = <GENOME>;
$line_before3 = <GENOME>;
$line=<GENOME>;
$line_after1=<GENOME>;
$line_after2=<GENOME>;
$line_after3=<GENOME>;
chomp($line_before1);
chomp($line_before2);
chomp($line_before3);
chomp($line);
chomp($line_after1);
chomp($line_after2);
chomp($line_after3);
my $line_local = uc($line_before1.$line_before2.$line_before3.$line.$line_after1.$line_after2.$line_after3);
#print "$line_local\n";
my @temp=split('',$line_local);
$position = $breakpoint3 % 50 - 1;
$coord = $position + length($line_before1) + length($line_before2) + length($line_before3);
$a = 0;
for($z=-140;$z<=140;$z++){
        $seq3[$a] = $temp[$coord+$z];
        $a++;
}
close GENOME;


if($chain3 eq '+'){
        printf "Right:%15s + 5' --- ",$geneRight;
        for($i=1;$i<140;$i++){
                print "$seq3[$i]";
        }
        for($i=140;$i<$#seq3;$i++){
                print  BOLD CYAN "$seq3[$i]", RESET;
		$chime .= $seq3[$i];
        }
        print " --- 3'\n";
}
elsif($chain3 eq '-'){
        printf "Right:%15s - 5' --- ",$geneRight;
        for($i=($#seq3-1);$i>140;$i--){
                print "$trans{$seq3[$i]}";
        }
        if($gjf eq "P9_BCL2"){
        	for($i=140;$i>=126;$i--){
			print BOLD CYAN "$trans{$seq3[$i]}", RESET;
			$chime .= $trans{$seq3[$i]};
		}
		for($i=125;$i>=123;$i--){
			print "-";
			$chime .= "-";
		}
		for($i=125;$i>=4;$i--){
			print BOLD CYAN "$trans{$seq3[$i]}", RESET;
			$chime .= $trans{$seq3[$i]};
		}
        }
        else{
		for($i=140;$i>=1;$i--){
			print BOLD CYAN "$trans{$seq3[$i]}", RESET;
			$chime .= $trans{$seq3[$i]};
		}
	}
        print " --- 3'\n";
}
else{
        print "Warning!!! Wrong chain signal\n";
}

printf "                               "; 
$fit = "";
for($i=0;$i<$insertion;$i++){
	$fit .= "-";
}
$chimeL = substr($chime,0,(140-$insertion-1));
$chimeR = substr($chime,(140-$insertion-1),140);
print BOLD YELLOW "$chimeL$fit$chimeR\n",RESET;
print "Split reads extracted from ";
print BOLD  RED "$sampleName[1]\n",RESET;
$seed1 = substr($chimeL,(105-$insertion),15);
$seed2 = substr($chimeL,(115-$insertion),15);
$seed3 = substr($chimeL,(125-$insertion),15);
$seed4 = substr($chimeR,0,15);
$seed5 = substr($chimeR,10,15);
$seed6 = substr($chimeR,20,15);
$blank = " " x 31;
$chime = $blank.$chime;


#@sorted_reads = sort{$a->[0] <=> $b->[0]} @records;

$correction_neg = 0;
$correction_pos = 0;
for($i=0;$i<=$#data;$i++){
	$match = index($data[$i][0],$seed1);
	if($match != -1){
		$data[$i][1] = 105 + 31 - $match -$insertion;
	#	print "seed0\n";
	}
}

for($i=0;$i<=$#data;$i++){
	$match = index($data[$i][0],$seed2);
	if($match != -1){
		$data[$i][1] = 115 + 31 - $match -$insertion;
	#	print "seed1\n";
	}
}

for($i=0;$i<=$#data;$i++){
        $match = index($data[$i][0],$seed3);
        if($match != -1){
                $data[$i][1] = 125 + 31 - $match -$insertion;
        #        print "seed2\n";
        }
}

for($i=0;$i<=$#data;$i++){
        $match = index($data[$i][0],$seed4);
        if($match != -1){
                $data[$i][1] = 140 + 31 - $match -1;
        #       print "seed3\n";
        }
}

for($i=0;$i<=$#data;$i++){
        $match = index($data[$i][0],$seed5);
        if($match != -1){
                $data[$i][1] = 150 + 31 - $match -1;
        #        print "seed4\n";
        }
}

for($i=0;$i<=$#data;$i++){
        $match = index($data[$i][0],$seed6);
        if($match != -1){
                $data[$i][1] = 160 + 31 - $match -1;
        #        print "seed4\n";
        }
}

for($i=0;$i<=$#data;$i++){
	my $start = " " x $data[$i][1];
        $data[$i][0] = $start.$data[$i][0];
}




my @zzz = grep { $_->[1] > -500 } @data;
@output = sort{$a->[1] <=> $b->[1]} @zzz;
@refer = split('',$chime);
for($i=0;$i<=$#output;$i++){
	@the_read = split('',$output[$i][0]);
	$xx = 170-$insertion;
	$lenEMP  = $output[$i][1];
	if($lenEMP  > $xx){
		$lenEMP  = $xx;
	}
	$empty = substr($output[$i][0],0,$lenEMP);
	print "$empty";
	#print "$output[$i][1]\n";
	for($j=$output[$i][1];$j<(170-$insertion);$j++){
		if($the_read[$j] ne $refer[$j]){
			print BOLD BLUE "$the_read[$j]",RESET;
		}
		else{
			print BOLD MAGENTA "$the_read[$j]", RESET;
		}
	}
	for($j=(170-$insertion);$j<170;$j++){
		print BOLD BLUE "$the_read[$j]",RESET;
	}
	for($j=170;$j<=$#the_read;$j++){
                if($the_read[$j] ne $refer[$j-$insertion]){
                        print BOLD BLUE "$the_read[$j]",RESET;
                }
		else{
			print BOLD CYAN "$the_read[$j]", RESET;
		}
	}
	print "\n";
}
print "\n";
