head -n1  header.txt > CNV_HG.seg
for i in {P2,P6,P7,P8};do
	cnvkit.py fix Coverage/$i\_HG.targetcoverage.cnn Coverage/$i\_HG.antitargetcoverage.cnn $i\_reference.cnn --no-edge -o subtraction/$i\_HG.cnr

	head -n1 subtraction/$i\_HG.cnr > test_HG.cnr
	tail -n+2 subtraction/$i\_HG.cnr > tmp_HG.bed
	bedtools subtract -a tmp_HG.bed -b ../Mask.bed >> test_HG.cnr

        cnvkit.py segment test_HG.cnr -p 20 --drop-low-coverage -m cbs -t 0.0000001 -o segmentation/$i\_HG.cbs.cns
        cnvkit.py call segmentation/$i\_HG.cbs.cns -o results/$i\_HG.cbs.call.cns
        cnvkit.py export seg results/$i\_HG.cbs.call.cns -o results/$i\_HG.cbs.seg

        tail -n+2 results/$i\_HG.cbs.seg >> CNV_HG.seg
done



for i in {P1,P3,P4,P5,P9,DHL38,DHL42};do
        cnvkit.py fix Coverage/$i\_HG.targetcoverage.cnn Coverage/$i\_HG.antitargetcoverage.cnn Super_normal_reference.cnn --no-edge -o subtraction/$i\_HG.cnr

        head -n1 subtraction/$i\_HG.cnr > test_HG.cnr
        tail -n+2 subtraction/$i\_HG.cnr > tmp_HG.bed
        bedtools subtract -a tmp_HG.bed -b ../Mask.bed >> test_HG.cnr

        cnvkit.py segment test_HG.cnr -p 20 --drop-low-coverage -m cbs -t 0.0000001 -o segmentation/$i\_HG.cbs.cns
        cnvkit.py call segmentation/$i\_HG.cbs.cns -o results/$i\_HG.cbs.call.cns
        cnvkit.py export seg results/$i\_HG.cbs.call.cns -o results/$i\_HG.cbs.seg

        tail -n+2 results/$i\_HG.cbs.seg >> CNV_HG.seg
done
