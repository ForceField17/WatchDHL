cat snv_FL_somatic_Mut.specific.txt           >  FL_somatic_Mut.specific.txt
tail -n+2 indel_FL_somatic_Mut.specific.txt  | cut -f 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113 >>  FL_somatic_Mut.specific.txt
cat snv_DHL_somatic_Mut.specific.txt          > DHL_somatic_Mut.specific.txt
tail -n+2 indel_DHL_somatic_Mut.specific.txt | cut -f 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113 >> DHL_somatic_Mut.specific.txt

cut -f 1,2,3,4,5,52,105 DHL_somatic_Mut.specific.txt > DHL.specific.bed

cat FL_somatic_Mut.specific.txt | awk -F"\t" '(NR==1 || $1!="P8" || ($105=="SNV" && $109>=4 && $50>=10) || ($105!="SNV" && $109>=5 && $109/($109+$108)*100>=10 ))'                   | cut -f 1,2,3,4,5,50,105 | sed s/P8/P8.a/g  > FL.specific.bed
cat FL_somatic_Mut.specific.txt | awk -F"\t" '( ($1=="P8" && $105=="SNV" && $51>=20 && $51!="NA" && $111>=5) || ($1=="P8" && $105!="SNV" && $111>=5 && $111/($110+$111)*100>=20) )'      | cut -f 1,2,3,4,5,51,105 | sed s/P8/P8.b/g >> FL.specific.bed


head -n1 DHL.specific.bed > DHL.specific.seg
head -n1 FL.specific.bed > FL.specific.seg
tail -n+2 DHL.specific.bed | awk -F"\t" '{if($7=="SNV"){print $1".DHL\t"$2"\t"$3-1"\t"$3"\t"$4"\t"$5"\t1"}else{print $1".DHL\t"$2"\t"$3"\t"$3+length($4)+length($5)-2"\t"$4"\t"$5"\t1"}}' | sort >> DHL.specific.seg
tail -n+2 FL.specific.bed  | awk -F"\t" '{if($7=="SNV"){print  $1".FL\t"$2"\t"$3-1"\t"$3"\t"$4"\t"$5"\t1"}else{print  $1".FL\t"$2"\t"$3"\t"$3+length($4)+length($5)-2"\t"$4"\t"$5"\t1"}}' | sort >> FL.specific.seg
