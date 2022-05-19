
for i in {P1,P2,P3,P4,P5,P6,P7,P9};do
   tail -n+2 ../$i/merge.report.coding.PDfilter.txt | awk -F"\t" '( $47<=1 )'  | awk -v name=$i  '{print name"\t"$0}' > temp_$i\.txt
   perl rmComma.pl temp_$i\.txt > tmp_$i\.txt
   perl rmPlus.pl tmp_$i\.txt  > ./results_Aug2021/$i\.raw.txt &
done

i=P8
   tail -n+2 ../$i/merge.report.coding.PDfilter.txt | awk -F"\t" '( $48<=1 )'  | awk -v name=$i  '{print name"\t"$0}' > temp_$i\.txt
   perl rmComma.pl temp_$i\.txt > tmp_$i\.txt
   perl rmPlus.pl tmp_$i\.txt  > ./results_Aug2021/$i\.raw.txt &
wait

for i in {P1,P2,P3,P4,P5,P6,P7,P8,P9};do
	rm tmp_$i\.txt &
	rm temp_$i\.txt &
done
