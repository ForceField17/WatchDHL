j=$1
tail -n+2 ../step4.FL.txt | awk -v CASE=$j '($53>=CASE || ($54>=CASE && $4=="P8") )' |  cut -f  1,2,4,7,8,108,122,123,128,129,133 | awk -F"\t" '{if($10~/SE/){print $0"\tSE"}else if($7~/H3K/){print $0"\tH3K27"}else{print $0"\tzNone"}}' | awk -F"\t" '{if($9~/Both/){print $0"\tBoth"}else if($7~/H3K/){print $0"\tH3K27"}else if($8~/H3K/){print $0"\tH3K4"}else{print $0"\tzNone"}}' > tmp1_FL
tail -n+2 tmp1_FL > tmp2_FL 
paste tmp1_FL tmp2_FL > tmp3_FL
cat tmp3_FL | awk '{if($1==$14){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$15-$2}}' > kateagis_FL.txt 
cat kateagis_FL.txt | grep -v chrY | sed s/chrX/23/g | sed s/chr//g | sort -nk1 -nk2 > sorted.kates.FL
perl ./absolute_pos.sh sorted.kates.FL > kateagis_FL_$j\.final.txt 
