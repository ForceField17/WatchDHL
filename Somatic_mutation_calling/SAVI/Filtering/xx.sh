
for i in {P3,P4,P5};do
   cat ../../../Patients_2022/$i\.raw.txt  | sort | uniq | awk -v name=$i  '{print name"\t"$0}' > temp_$i\.txt
   perl rmComma.pl temp_$i\.txt > tmp_$i\.txt
   perl rmPlus.pl tmp_$i\.txt  > ./$i\.raw.txt &
done

for i in {P3,P4,P5};do
	rm tmp_$i\.txt &
	rm temp_$i\.txt &
done
