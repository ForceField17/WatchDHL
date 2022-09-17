i=$1
j=$2
tail -n+2 ../../step4.$i\.txt | awk -v CASE=$j '($4==CASE)' |  cut -f  1,2,4,7,8,108,122,123,128,129,133 | awk -F"\t" '{if($10~/SE/){print $0"\tSE"}else if($7~/H3K/){print $0"\tH3K27"}else{print $0"\tzNone"}}' | awk -F"\t" '{if($9~/Both/){print $0"\tBoth"}else if($7~/H3K/){print $0"\tH3K27"}else if($8~/H3K/){print $0"\tH3K4"}else{print $0"\tzNone"}}' > tmp1_$i
tail -n+2 tmp1_$i > tmp2_$i 
paste tmp1_$i tmp2_$i > tmp3_$i
cat tmp3_$i | awk '{if($1==$14){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$15-$2}}' > kateagis_$i\.txt 
cat kateagis_$i\.txt | grep -v chrY | sed s/chrX/23/g | sed s/chr//g | sort -nk1 -nk2 > sorted.kates.$i
perl ./absolute_pos.sh sorted.kates.$i > kateagis_$i\_$j\.final.txt 
