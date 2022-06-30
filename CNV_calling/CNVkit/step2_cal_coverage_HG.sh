for i in {P1,P2,P3,P4,P5,P6,P7,P8,P9};do
	cnvkit.py coverage /scratch/PI/jgwang/dsongad/DLBCL/raw_data/results_T_DNA/$i\_HG.sorted.MD.bam DHL_target.bed     -p 20 -o Coverage/$i\_HG.targetcoverage.cnn &
	cnvkit.py coverage /scratch/PI/jgwang/dsongad/DLBCL/raw_data/results_T_DNA/$i\_HG.sorted.MD.bam DHL_antitarget.bed -p 2  -o Coverage/$i\_HG.antitargetcoverage.cnn &
wait
done
wait

for i in {DHL38,DHL42};do
        cnvkit.py coverage /scratch/PI/jgwang/dsongad/DLBCL/raw_data/results_T_DNA/$i\.sorted.MD.bam DHL_target.bed     -p 40 -o Coverage/$i\_HG.targetcoverage.cnn &
        cnvkit.py coverage /scratch/PI/jgwang/dsongad/DLBCL/raw_data/results_T_DNA/$i\.sorted.MD.bam DHL_antitarget.bed -p 2  -o Coverage/$i\_HG.antitargetcoverage.cnn 
done
wait
