cat ../filtered.P345.txt | awk -F"\t" '(NR==1 || ($106==0 && $50==0 && $110>=3 && $52>=15 && $46>=12))' | awk -F"\t" '(NR==1 || (($109+1)/($110+1))<2 )' | awk -F"\t" '(NR==1 || ($52>=20 && $110>=5) || $1~/^P[3]/ )'   > DHL_P345_Mut.specific.txt
cat ../filtered.P345.txt | awk -F"\t" '(NR==1 || ($106>=3 && $50>=15 && $110==0 && $52==0 && ($48>=12 || ($48>=7 && $1=="P3") )))' | awk -F"\t" '(NR==1 || (($105+1)/($106+1))<2 )' | awk -F"\t" '(NR==1 || ($50>=20 && $106>=5) || $1~/^P[3]/ )'   > FL_P345_Mut.specific.txt

cut -f 1,2,3,4,5,52,100 DHL_P345_Mut.specific.txt > DHL.P345.bed

cut -f 1,2,3,4,5,50,100 FL_P345_Mut.specific.txt  > FL.P345.bed

cat DHL.P345.bed | awk -F"\t" '{if($7=="SNP"){print $1".DHL\t"$2"\t"$3-1"\t"$3"\t"$4"\t"$5"\t"$6}else{print $1".DHL\t"$2"\t"$3"\t"$3+length($4)+length($5)-2"\t"$4"\t"$5"\t"$6}}' > DHL.P345.seg
cat FL.P345.bed  | awk -F"\t" '{if($7=="SNP"){print $1".FL\t"$2"\t"$3-1"\t"$3"\t"$4"\t"$5"\t"$6}else{print $1".FL\t"$2"\t"$3"\t"$3+length($4)+length($5)-2"\t"$4"\t"$5"\t"$6}}' > FL.P345.seg
