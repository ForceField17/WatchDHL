for i in {P1,P2,P4,P5,P6,P7,P8,P9,P8_2};do
	cnvkit.py coverage /scratch/PI/jgwang/dsongad/DLBCL/raw_data/results_T_DNA/$i\_LG.sorted.MD.bam DHL_target.bed     -p 20 -o Coverage/$i\_LG.targetcoverage.cnn &
	cnvkit.py coverage /scratch/PI/jgwang/dsongad/DLBCL/raw_data/results_T_DNA/$i\_LG.sorted.MD.bam DHL_antitarget.bed -p 2  -o Coverage/$i\_LG.antitargetcoverage.cnn &
wait
done
wait
