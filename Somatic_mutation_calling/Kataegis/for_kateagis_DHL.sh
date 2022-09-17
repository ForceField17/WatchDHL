j=$1
tail -n+2 ../step4.DHL.txt | awk -v CASE=$j '($55>=CASE )' |  cut -f  1,2,4,7,8,108,122,123,128,129,133 | awk -F"\t" '{if($10~/SE/){print $0"\tSE"}else if($7~/H3K/){print $0"\tH3K27"}else{print $0"\tzNone"}}' | awk -F"\t" '{if($9~/Both/){print $0"\tBoth"}else if($7~/H3K/){print $0"\tH3K27"}else if($8~/H3K/){print $0"\tH3K4"}else{print $0"\tzNone"}}' > tmp1_DHL
tail -n+2 tmp1_DHL > tmp2_DHL 
paste tmp1_DHL tmp2_DHL > tmp3_DHL
cat tmp3_DHL | awk '{if($1==$14){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$15-$2}}' > kateagis_DHL.txt 
cat kateagis_DHL.txt | grep -v chrY | sed s/chrX/23/g | sed s/chr//g | sort -nk1 -nk2 > sorted.kates.DHL
perl ./absolute_pos.sh sorted.kates.DHL > kateagis_DHL_$j\.final.txt 
