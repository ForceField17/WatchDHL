#mkdir heatmap
i=P8
head -n1  ../../step4.FL.csv | cut -f 4,5,6,7,8,47,48,49,50,52,53,54,55,75,84,93,108 | awk '{print $0"\tgermline\tFL.a\tFL.b\tDHL"}' | sed s/\#//g > $i.mutation.heatmap.txt 
cat ../../step4.*csv | sort | uniq | awk -F"\t" -v xxx=$i '($4==xxx && $47>=15 && $48>=15 && $49>=15 && $50>=15)' | cut -f 4,5,6,7,8,47,48,49,50,52,53,54,55,75,84,93,108 | awk '{if($10>=10){print $0"\t1"}else{print $0"\t0"}}'   | awk '{if($11>=2 && $14>=2){print $0"\t1"}else{print $0"\t0"}}'  | awk '{if($12>=2 && $15>=2 && ($12>=20 || $13>=20)){print $0"\t1"}else{print $0"\t0"}}'  | awk '{if($13>=2 && $16>=2 && ($12>=20 || $13>=20)){print $0"\t1"}else{print $0"\t0"}}'  >> $i.mutation.heatmap.txt


#P1
for i in {P6,P7};do
head -n1  ../../step4.FL.csv | cut -f 4,5,6,7,8,47,48,50,52,53,55,108 | awk '{print $0"\tgermline\tFL\tDHL"}' | sed s/\#//g > $i.mutation.heatmap.txt
cat ../../step4.*csv | sort | uniq | awk -F"\t" -v xxx=$i '($4==xxx && $47>=15 && $48>=15 && $50>=15)' | cut -f 4,5,6,7,8,47,48,50,52,53,55,108 | awk '{if($9>=10){print $0"\t1"}else{print $0"\t0"}}'   | awk '{if($10>=2){print $0"\t1"}else{print $0"\t0"}}'  | awk '{if($11>=2){print $0"\t1"}else{print $0"\t0"}}' >> $i.mutation.heatmap.txt
done

i=P9
head -n1  ../../step4.FL.csv | cut -f 4,5,6,7,8,47,48,50,52,53,55,108 | awk '{print $0"\tgermline\tFL\tDHL"}' | sed s/\#//g > $i.mutation.heatmap.txt
cat ../../step4.*csv | sort | uniq | awk -F"\t" -v xxx=$i '($4==xxx && $47>=15 && $48>=15 && $50>=15)' | cut -f 4,5,6,7,8,47,48,50,52,53,55,108 | awk '{if($9>=10){print $0"\t1"}else{print $0"\t0"}}'   | awk '{if($10>=2){print $0"\t1"}else{print $0"\t0"}}'  | awk '{if($11>=2){print $0"\t1"}else{print $0"\t0"}}' >> $i.mutation.heatmap.txt


i=P1
head -n1  ../../step4.FL.csv | cut -f 4,5,6,7,8,47,48,50,52,53,55,108 | awk '{print $0"\tgermline\tFL\tDHL"}' | sed s/\#//g > $i.mutation.heatmap.txt
cat ../../step4.*csv | sort | uniq | awk -F"\t" -v xxx=$i '($4==xxx && $47>=15 && $48>=15 && $50>=15)' | cut -f 4,5,6,7,8,47,48,50,52,53,55,108 | awk '{if($9>=10){print $0"\t1"}else{print $0"\t0"}}'   | awk '{if($10>=2){print $0"\t1"}else{print $0"\t0"}}'  | awk '{if($11>=2){print $0"\t1"}else{print $0"\t0"}}' >> $i.mutation.heatmap.txt

i=P2
head -n1  ../../step4.FL.csv | cut -f 4,5,6,7,8,47,48,50,52,53,55,108 | awk '{print $0"\tgermline\tFL\tDHL"}' | sed s/\#//g > $i.mutation.heatmap.txt
cat ../../step4.*csv | sort | uniq | awk -F"\t" -v xxx=$i '($4==xxx && $47>=15 && $48>=15 && $50>=15)' | cut -f 4,5,6,7,8,47,48,50,52,53,55,108 | awk '{if($9>=10){print $0"\t1"}else{print $0"\t0"}}'   | awk '{if($10>=2){print $0"\t1"}else{print $0"\t0"}}'  | awk '{if($11>=2){print $0"\t1"}else{print $0"\t0"}}' >> $i.mutation.heatmap.txt
