cat kateagis_FL_20.final.txt | awk '(NR==1 || $14<=5000)' | awk '{print $3"\tchr"$1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t1"}' | grep -v chr23 > FL_20_mutation_clusters.seg
