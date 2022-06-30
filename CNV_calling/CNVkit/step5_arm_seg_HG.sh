head -n1  header.txt > Arm_HG.seg
for i in {P1,P2,P3,P4,P5,P6,P7,P8,P9,DHL38,DHL42};do

	head -n1 subtraction/$i\_HG.cnr > test_HG.cnr
	tail -n+2 subtraction/$i\_HG.cnr > tmp_HG.bed
	bedtools subtract -a tmp_HG.bed -b ../Mask.bed >> test_HG.cnr

        cnvkit.py segment test_HG.cnr -p 20 --drop-low-coverage -m none -t 0.0000001 -o segmentation/$i\_HG.Arm.cns
        cnvkit.py call segmentation/$i\_HG.Arm.cns -o results/$i\_HG.Arm.call.cns
        cnvkit.py export seg results/$i\_HG.Arm.call.cns -o results/$i\_HG.Arm.seg

        tail -n+2 results/$i\_HG.Arm.seg >> Arm_HG.seg
done


