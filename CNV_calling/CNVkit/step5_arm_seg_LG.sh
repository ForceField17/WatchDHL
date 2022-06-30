head -n1  header.txt > Arm_LG.seg
for i in {P1,P2,P3,P4,P5,P6,P7,P8,P8_2,P9};do

	head -n1 subtraction/$i\_LG.cnr > test_LG.cnr
	tail -n+2 subtraction/$i\_LG.cnr > tmp_LG.bed
	bedtools subtract -a tmp_LG.bed -b ../Mask.bed >> test_LG.cnr

        cnvkit.py segment test_LG.cnr -p 20 --drop-low-coverage -m none -t 0.000001 -o segmentation/$i\_LG.Arm.cns
        cnvkit.py call segmentation/$i\_LG.Arm.cns -o results/$i\_LG.Arm.call.cns
        cnvkit.py export seg results/$i\_LG.Arm.call.cns -o results/$i\_LG.Arm.seg

        tail -n+2 results/$i\_LG.Arm.seg >> Arm_LG.seg
done


