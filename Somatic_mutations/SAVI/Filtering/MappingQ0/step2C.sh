
for i in {P1,P2,P3,P4,P5,P6,P7,P8,P9};do
	sh MQ0B.sh $i\.2A.txt $i\.2B.txt $i\_LG.sorted.MD.bam &
done
wait

for i in {P1,P2,P3,P4,P5,P6,P7,P8,P9};do
        sh MQ0B.sh $i\.2B.txt $i\.2C.txt $i\_HG.sorted.MD.bam &
done
wait

