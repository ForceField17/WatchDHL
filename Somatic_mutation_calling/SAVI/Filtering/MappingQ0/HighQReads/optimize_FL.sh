cat filtered.all_somatic.txt | awk '($1=="P2" && ($106>=2 || $105>=2) )' | cut -f 1,2,3,4,5 | awk '{print $1"\t"$2"\t"$3"\t"$3+1"\t"$4"\t"$5"\t1"}' > P2.FL.seg


cat filtered.all_somatic.txt | awk '($1~/^P[289]/ && ($106+$105>=4) )' | cut -f 1,2,3,4,5 | awk '{print $1"\t"$2"\t"$3"\t"$3+1"\t"$4"\t"$5"\t1"}' > P289.FL.seg
