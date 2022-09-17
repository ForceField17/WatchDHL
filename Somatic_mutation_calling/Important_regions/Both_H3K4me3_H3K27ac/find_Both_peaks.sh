cut -f 4 ../../../CellLineEpi/Both.H3K4me3_and_H3K27ac.bed  | sort > tmp
Total_FL=`tail -n+2 ../../step4.FL.csv | wc -l`
Total_DHL=`tail -n+2 ../../step4.DHL.csv | wc -l`
cut -f 128 ../../step4.DHL.csv | tail -n+2 | grep -v "\." | sort | uniq -c | awk -v xxx=$Total_DHL '{print $2"\t"($1/xxx)*100"\t"$1}' | sort > DHL_SE.txt
cut -f 128 ../../step4.FL.csv  | tail -n+2 | grep -v "\." | sort | uniq -c | awk -v xxx=$Total_FL '{print $2"\t"($1/xxx)*100"\t"$1}' | sort > FL_SE.txt
join tmp FL_SE.txt -a1 | awk '{if($2==""){print $1"\t0\t0"}else{print $1"\t"$2"\t"$3}}' | sort > tmpFL
join tmpFL DHL_SE.txt -a1 | awk '{if($4==""){print $1"\t"$2"\t"$3"\t0\t0"}else{print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' | sort > tmpDHL
cat tmpDHL >tmp2 #| awk '($3>0 || $5>0)' > tmp2
sort-bed ../../counting_mutations_in_BINs/summary_in_gene_level/gene.bed | cut -f 1,2,3,6 > hg38.gene.bed
bedtools closest -a ../../../CellLineEpi/Both.H3K4me3_and_H3K27ac.bed  -b hg38.gene.bed -t first -d | awk '{print $4"\t"$3-$2"\t"$8}' | sort > List.txt
join tmp2 List.txt | sed s/" "/\\t/g | awk '{print $0"\t"$2*1000000/$6"\t"$4*1000000/$6}' | awk '{print $0"\t"$4-$2}' > MutationPercent.txt

cat MutationPercent.txt | awk '{print $0"\t"$5-$3}' | sort -nrk4 | awk '{count++}{print $0"\t"count}' | sort -nrk5 | awk '{count++}{print $0"\t"count}' | sort -nrk11 | awk '{count++}{print $0"\t"count}' > results.txt
