mkdir Coverage
mkdir segmentation
mkdir results
mkdir GISTIC
mkdir subtraction
cnvkit.py reference   -f /scratch/PI/jgwang/dsongad/software/my_lib/star_genome_d1_vd1_gtfv22/GRCh38.d1.vd1.fa --no-edge -o flat_reference.cnn -t DHL_target.bed -a DHL_antitarget.bed &
for i in {P1,P2,P6,P7,P8,P9,Super_normal};do	
	cnvkit.py coverage /scratch/PI/jgwang/dsongad/DLBCL/raw_data/results_B_DNA/$i\.bam DHL_target.bed -p 20 -o $i\.targetcoverage.cnn &
	cnvkit.py coverage /scratch/PI/jgwang/dsongad/DLBCL/raw_data/results_B_DNA/$i\.bam DHL_antitarget.bed -p 20 -o $i\.antitargetcoverage.cnn &
	wait
	cnvkit.py reference   -f /scratch/PI/jgwang/dsongad/software/my_lib/star_genome_d1_vd1_gtfv22/GRCh38.d1.vd1.fa --no-edge -o $i\_reference.cnn $i\.*coverage.cnn &
done
wait
