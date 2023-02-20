#!/usr/bin/perl

#please run this job in a new clear fold
die "usage: perl checkReads.mpi.pl <preprocessed_savi_table> <TumorBam_Folder> \n" if @ARGV!= 2;

my %trans;


$target = $ARGV[0];
$foldT = $ARGV[1];
open FILE, $target;
#$head = <FILE>;
chomp($head);
#print "$head\tref\talt\tseqL\tseqR\tNgoodRefFather\tNgoodAltFather\tAFgoodFather\tNgoodRefMother\tNgoodAltMother\tAFgoodMother\tNgoodRefBlood\tNgoodAltBlood\tAFgoodBlood\tNgoodRefTumor\tNgoodAltTumor\tAFgoodTumor\n";
$i=0;
while($line = <FILE>){
	chomp($line);
	$content[$i] = $line;
	@temp=split('\t',$line);
	$chr = $temp[1];
	$pos = $temp[2];
	$ID = $temp[0];
	$obj1 = $pos-20;
	$obj2 = $pos+20;	
	my $left = `samtools faidx /scratch/PI/jgwang/dsongad/software/my_lib/star_genome_d1_vd1_gtfv22/GRCh38.d1.vd1.fa $chr\:$obj1\-$pos | tail -n +2`;
	$left =~ s/\n//g;
	chomp($left);
        my $right = `samtools faidx /scratch/PI/jgwang/dsongad/software/my_lib/star_genome_d1_vd1_gtfv22/GRCh38.d1.vd1.fa $chr\:$pos\-$obj2 | tail -n +2`;
        $right =~ s/\n//g;
        chomp($right);
	$REF = $temp[3];
	$ALT = $temp[4];
	$seqL[$i] = $left;
	$seqR[$i] = $right;
	$fileFLa = $temp[101];
	$fileFLb = $temp[102];
	$fileDHL = $temp[103];	
	$mutType = $temp[104];
	#print "$chr\t$pos\t$ID\t$REF\t$ALT\t$fold\t$fileFLb\t$mutType\n";
    if($mutType eq 'Ins'){	
		$len = length($ALT)-1;

            $mapping1_2 = `samtools view -b /scratch/PI/jgwang/dsongad/DLBCL/raw_data/results_B_DNA/Super_normal.bam $chr\:$pos\-$pos | samtools calmd -e - /scratch/PI/jgwang/dsongad/software/my_lib/star_genome_d1_vd1_gtfv22/GRCh38.d1.vd1.fa | grep -v "^@" | awk -F"\t" '(\$5>=30)' | awk -v pos=$pos 'BEGIN {OFS = FS = "\t" } ; {theX=(pos-\$4+1); theY=$len; if(\$6 !~ theX"M" ) print \$4 }' | wc -l`;
            $mapping1_2 =~ s/\n//g;
            $Nref_2[$i] = $mapping1_2;
            $mapping2_2 = `samtools view -b /scratch/PI/jgwang/dsongad/DLBCL/raw_data/results_B_DNA/Super_normal.bam $chr\:$pos\-$pos | samtools calmd -e - /scratch/PI/jgwang/dsongad/software/my_lib/star_genome_d1_vd1_gtfv22/GRCh38.d1.vd1.fa | grep -v "^@" | awk -F"\t" '(\$5>=30)' |  awk -v pos=$pos 'BEGIN {OFS = FS = "\t" } ; {theX=(pos-\$4+1); theY=$len; if(\$6 ~ theX"M"theY"I" ) print \$4 }' | wc -l`;
            $mapping2_2 =~ s/\n//g;
            $Nalt_2[$i] = $mapping2_2

	} 
    elsif($mutType eq 'Del'){
        $len = length($REF)-1;


            $mapping1_2 = `samtools view -b /scratch/PI/jgwang/dsongad/DLBCL/raw_data/results_B_DNA/Super_normal.bam $chr\:$pos\-$pos | samtools calmd -e - /scratch/PI/jgwang/dsongad/software/my_lib/star_genome_d1_vd1_gtfv22/GRCh38.d1.vd1.fa | grep -v "^@" | awk -F"\t" '(\$5>=30)' | awk -v pos=$pos 'BEGIN {OFS = FS = "\t" } ; {theX=(pos-\$4+1); theY=$len; if(\$6 !~ theX"M" ) print \$4 }' | wc -l`;
            $mapping1_2 =~ s/\n//g;
            $Nref_2[$i] = $mapping1_2;
            $mapping2_2 = `samtools view -b /scratch/PI/jgwang/dsongad/DLBCL/raw_data/results_B_DNA/Super_normal.bam $chr\:$pos\-$pos | samtools calmd -e - /scratch/PI/jgwang/dsongad/software/my_lib/star_genome_d1_vd1_gtfv22/GRCh38.d1.vd1.fa | grep -v "^@" | awk -F"\t" '(\$5>=30)' |  awk -v pos=$pos -v CHR=$chr 'BEGIN {OFS = FS = "\t" } ; {theX=(pos-\$4+1); theY=$len; POS=pos+theY+1; if(\$6 ~ theX"M"theY"D" || (\$6 ~ theX"M[0-9]+S" && \$12 ~ "SA:Z:"CHR","POS"," ) ) print \$4 }'  | wc -l`;
            $mapping2_2 =~ s/\n//g;
            $Nalt_2[$i] = $mapping2_2;

    } 
	else{

        $Nref_2[$i] = "-";
        $Nalt_2[$i] = "-";

	}


	print "$content[$i]\t$Nref_2[$i]\t$Nalt_2[$i]\n";
	$i++;
}
close FILE;

