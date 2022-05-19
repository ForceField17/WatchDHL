cat ../filtered.all_somatic.txt | awk -F"\t" '(NR==1 || ($110>=4 && $52>=5))' | awk -F"\t" '(NR==1 || (($109+1)/($110+1))<2 )' | awk -F"\t" '(NR==1 || ($52>=20 && $110>=5) || $1~/^P[289]/ )'   > DHL_somatic_Mut.specific.txt
cat ../filtered.all_somatic.txt | awk -F"\t" '(NR==1 || ($106>=4 && $50>=5)  || ($51>=20 && $108>=5 && $108!="-") )' | awk -F"\t" '(NR==1 || (($105+1)/($106+1))<2  || ( (($107+1)/($108+1))<2  && $51>=20 && $108>=5 && $108!="-")  )' | awk -F"\t" '(NR==1 || ($50>=20 && $106>=5) || ($51>=20 && $108>=5 && $108!="-") || $1~/^P[289]/  )' > FL_somatic_Mut.specific.txt

cut -f 1,2,3,4,5,52,100 DHL_somatic_Mut.specific.txt > DHL.specific.bed
cut -f 1,2,3,4,5,52,100 SV.txt | awk '($6>=10 && $6!="NA")' >> DHL.specific.bed

cat FL_somatic_Mut.specific.txt | awk -F"\t" '(NR==1 || ($106>=4 && $50>=5))' | cut -f 1,2,3,4,5,50,100  > FL.specific.bed
cut -f 1,2,3,4,5,51,100 FL_somatic_Mut.specific.txt | awk '($6>=20 && $6!="NA")' | grep -v CaseID |  sed s/P8/P8.b/g >> FL.specific.bed

cut -f 1,2,3,4,5,50,100 SV.txt | awk '($6>=10 && $6!="NA")' >> FL.specific.bed
cut -f 1,2,3,4,5,51,100 SV.txt | awk '($6>=10 && $6!="NA")' |  sed s/P8/P8.b/g  >> FL.specific.bed

cat DHL.specific.bed | awk -F"\t" '{if($7=="SNP"){print $1".DHL\t"$2"\t"$3-1"\t"$3"\t"$4"\t"$5"\t"$6}else{print $1".DHL\t"$2"\t"$3"\t"$3+length($4)+length($5)-2"\t"$4"\t"$5"\t"$6}}' > DHL.specific.seg
cat FL.specific.bed | awk -F"\t" '{if($7=="SNP"){print $1".FL\t"$2"\t"$3-1"\t"$3"\t"$4"\t"$5"\t"$6}else{print $1".FL\t"$2"\t"$3"\t"$3+length($4)+length($5)-2"\t"$4"\t"$5"\t"$6}}' > FL.specific.seg
