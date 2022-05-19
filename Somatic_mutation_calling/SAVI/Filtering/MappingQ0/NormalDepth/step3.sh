
for i in {P2,P6};do
	sh MQ0B.sh $i\.2C.txt $i\.3.txt $i\.bam &
done
wait

