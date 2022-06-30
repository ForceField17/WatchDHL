#

# Author: @naveedishaque
# Date 19 Feb 2019

use strict;
use warnings;

# read file or exit
my $file = shift or die("ERROR: expected 'scores.gistic' input file\n\n");

my $del_genes_file = shift or die("ERROR: expected 'del_genes.conf20.txt' input file\n\n");
die "ERROR: cannot open file '$file'\n" unless -e $file;

my $amp_genes_file = shift or die("ERROR: expected 'amp_genes.conf20.txt' input file\n\n");
die "ERROR: cannot open file '$file'\n" unless -e $file;

my %offsets;

# open file and parse header. Exit if header is unexpected

open (my $gistic_fh, "<$file");
my @gistic_array = <$gistic_fh>;
chomp(@gistic_array);
close($gistic_fh);

my $expected_header="Type\tChromosome\tStart\tEnd\t-log10(q-value)\tG-score\taverage amplitude\tfrequency";

if  ($gistic_array[0] ne $expected_header){
  die "ERROR: input file '$file' does not have the expected header....\n\n\tObserved header: '$gistic_array[0]'\n\tExpected: '$expected_header'\n\n";
}

my $num_results = scalar(@gistic_array);

# iterate over the file, and produce a simple datafame that can be plotted as a polygon using R
# save Chr end possitions for plotting!

my @chr_lines;

my ($last_chr, $last_end, $offset, $last_type) = (" 1",0,0,"Amp");
$offsets{$last_chr}=$offset;

print "Type\tChromosome\tPosition\tPosition_offset\t-log10(q-value)\tG-score\taverage_amplitude\tfrequency\n";

foreach my $gistic_line(@gistic_array){
  my @l = split ("\t",$gistic_line);
  next if $l[0] eq "Type";
  if ($l[4] > 40){
    warn "WARNING: encountered a region with suspiciously high signficance... possiblly a fasle positive. Skipping line:\n\t$gistic_line\n\n";
	next;
  }

  if ($l[1] ne $last_chr){
    push(@chr_lines,"Chr\t$last_chr\t$last_end\t".($last_end+$offset)."\t".($offset+$last_end/2)."\tNA\tNA\tNA\n") if ($last_type eq "Amp");
    $offsets{$last_chr}=$offset;
    $offset = $offset + $last_end;
    $offset = 0 if ($l[0] ne $last_type);
	$last_end = 0;
  }
  
  # set regions to 0 if they are skipped by gistic
  if($l[2] > ($last_end + 1)){
    my $zero_start =  $last_end + 1;
    my $zero_end = $l[2] -1;
    print "$l[0]\t$l[1]\t$zero_start\t".($zero_start+$offset)."\t$l[4]\t0\t0\t0\n";
    print "$l[0]\t$l[1]\t$zero_end\t".($zero_end+$offset)."\t$l[4]\t0\t0\t0\n";
  }

  print "$l[0]\t$l[1]\t$l[2]\t".($l[2]+$offset)."\t$l[4]\t$l[5]\t$l[6]\t$l[7]\n";
  print "$l[0]\t$l[1]\t$l[3]\t".($l[3]+$offset)."\t$l[4]\t$l[5]\t$l[6]\t$l[7]\n";

  ($last_chr, $last_end, $last_type) = ($l[1], $l[3], $l[0]);

}

print @chr_lines;

####

## parse gene lists

open(IN, "awk -F'\t' '{for (i=1; i<=NF; i++) {a[NR,i] = \$i}} NF>p { p = NF } END { for(j=1; j<=p; j++) { str=a[1,j]; for(i=2; i<=NR; i++){ str=str\"\t\"a[i,j]; } print str } }' $del_genes_file | cut -f 4- | grep chr | sed 's|chr||' | sed 's|-|:|' | sed 's|\t|:|' | sed 's|\t\t||g' | sed 's|\t\$||g' | sed 's|\t|,|g' | ") or die;
while (<IN>){
  chomp;
  my @l = split(":");
  $l[0] = " ".$l[0] if ($l[0] < 10);
  print "Gene_del\t$l[0]\t$l[1]\t".($l[1]+$offsets{$l[0]})."\t".$offsets{$l[0]}."\t0\t$l[3]\t0\n";
}
close(IN);

open(IN, "awk -F'\t' '{for (i=1; i<=NF; i++) {a[NR,i] = \$i}} NF>p { p = NF } END { for(j=1; j<=p; j++) { str=a[1,j]; for(i=2; i<=NR; i++){ str=str\"\t\"a[i,j]; } print str } }' $amp_genes_file | cut -f 4- | grep chr | sed 's|chr||' | sed 's|-|:|' | sed 's|\t|:|' | sed 's|\t\t||g'  | sed 's|\t\$||g' | sed 's|\t|,|g' | ") or die;
while (<IN>){
  chomp;
  my @l = split(":");
  $l[0] = " ".$l[0] if ($l[0] < 10);
  print "Gene_amp\t$l[0]\t$l[1]\t".($l[1]+$offsets{$l[0]})."\t".$offsets{$l[0]}."\t0\t$l[3]\t0\n";
}
close(IN);
