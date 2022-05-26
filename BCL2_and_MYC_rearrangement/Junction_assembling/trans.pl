#!/usr/bin/perl
#programmed by SONG,DONG 2019.12
#please run this job in a new clean fold
die "usage: perl fusion_visual.pl <report_file_1_line>\n" if @ARGV!= 1;
#use warnings;
$gjf = $ARGV[0];
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

$i = 0;
$j = 0;

@negChain5 = ();
@posChain5 = ();

chomp($gjf);
my @temp=split('',$gjf);
for($i=0;$i<=$#temp;$i++){
	print "$temp[$i]";
}
print "\n";
for($i=$#temp;$i>=0;$i--){
        print "$trans{$temp[$i]}";
}
print "\n";
