head -n1  header.txt > CNV_LG.seg
for i in {P2,P6,P7,P8,P8_2};do
	cnvkit.py fix Coverage/$i\_LG.targetcoverage.cnn Coverage/$i\_LG.antitargetcoverage.cnn $i\_reference.cnn --no-edge -o subtraction/$i\_LG.cnr

	head -n1 subtraction/$i\_LG.cnr > test_LG.cnr
	tail -n+2 subtraction/$i\_LG.cnr > tmp_LG.bed
	bedtools subtract -a tmp_LG.bed -b ../Mask.bed >> test_LG.cnr

        cnvkit.py segment test_LG.cnr -p 20 --drop-low-coverage -m cbs -t 0.000001 -o segmentation/$i\_LG.cbs.cns
        cnvkit.py call segmentation/$i\_LG.cbs.cns -o results/$i\_LG.cbs.call.cns
        cnvkit.py export seg results/$i\_LG.cbs.call.cns -o results/$i\_LG.cbs.seg

        tail -n+2 results/$i\_LG.cbs.seg >> CNV_LG.seg
done



for i in {P1,P4,P5,P9};do
        cnvkit.py fix Coverage/$i\_LG.targetcoverage.cnn Coverage/$i\_LG.antitargetcoverage.cnn Super_normal_reference.cnn --no-edge -o subtraction/$i\_LG.cnr

        head -n1 subtraction/$i\_LG.cnr > test_LG.cnr
        tail -n+2 subtraction/$i\_LG.cnr > tmp_LG.bed
        bedtools subtract -a tmp_LG.bed -b ../Mask.bed >> test_LG.cnr

        cnvkit.py segment test_LG.cnr -p 20 --drop-low-coverage -m cbs -t 0.000001 -o segmentation/$i\_LG.cbs.cns
        cnvkit.py call segmentation/$i\_LG.cbs.cns -o results/$i\_LG.cbs.call.cns
        cnvkit.py export seg results/$i\_LG.cbs.call.cns -o results/$i\_LG.cbs.seg

        tail -n+2 results/$i\_LG.cbs.seg >> CNV_LG.seg
done
