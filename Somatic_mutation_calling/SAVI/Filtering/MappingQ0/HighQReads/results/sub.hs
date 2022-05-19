cat ../filtered.all_somatic.txt | awk -F"\t" '(NR==1 || ($52>=20 && $110>=5))' > DHL_somatic_Mut.20percent_5reads.txt
cat ../filtered.all_somatic.txt | awk -F"\t" '(NR==1 || ($50>=20 && $106>=5) || ($51>=20 && $108>=5 && $108!="-") )' > FL_somatic_Mut.20percent_5reads.txt

cut -f 1,2,3,4,5,52,100 DHL_somatic_Mut.20percent_5reads.txt > DHL.20percent.bed
cut -f 1,2,3,4,5,52,100 SV.txt | awk '($6>=10 && $6!="NA")' >> DHL.20percent.bed

cut -f 1,2,3,4,5,50,100 FL_somatic_Mut.20percent_5reads.txt | awk '($6>=20 && $6!="NA")' > FL.20percent.bed
cut -f 1,2,3,4,5,51,100 FL_somatic_Mut.20percent_5reads.txt | awk '($6>=20 && $6!="NA")' |grep -v CaseID |  sed s/P8/P8.2/g >> FL.20percent.bed

cut -f 1,2,3,4,5,50,100 SV.txt | awk '($6>=10 && $6!="NA")' >> FL.20percent.bed
cut -f 1,2,3,4,5,51,100 SV.txt | awk '($6>=10 && $6!="NA")' |  sed s/P8/P8.2/g  >> FL.20percent.bed

cat DHL.20percent.bed | awk -F"\t" '{if($7=="SNP"){print $1".DHL\t"$2"\t"$3-1"\t"$3"\t"$4"\t"$5"\t"$6}else{print $1".DHL\t"$2"\t"$3"\t"$3+length($4)+length($5)-2"\t"$4"\t"$5"\t"$6}}' > DHL.20percent.seg
cat FL.20percent.bed | awk -F"\t" '{if($7=="SNP"){print $1".FL\t"$2"\t"$3-1"\t"$3"\t"$4"\t"$5"\t"$6}else{print $1".FL\t"$2"\t"$3"\t"$3+length($4)+length($5)-2"\t"$4"\t"$5"\t"$6}}' > FL.20percent.seg
