
for i in {P3,P4,P5};do
	sh MQ0B.sh $i\.2A.txt $i\.2B.txt $i\_LG.sorted.MD.bam &
done
wait

for i in {P3,P4,P5};do
        sh MQ0B.sh $i\.2B.txt $i\.2C.txt $i\_HG.sorted.MD.bam &
done
wait

