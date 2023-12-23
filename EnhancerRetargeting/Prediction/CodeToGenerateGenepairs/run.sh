tail -n+3 ./ENCFF156ECM.TAD.bed > TAD.bed
sort-bed TAD.bed >tmp
bedtools merge -i tmp | awk '{count+=1}{print $0"\tTAD_"count}' > TAD.sorted.bed

grep -v chrY TSS_hg38.bed | cut -f 6 | sort | uniq -c | awk '($1<2)' | awk '{print $2}' | sort > legal_genes.txt

cat TSS_hg38.bed | awk '{print $6"\t"$0}' | sort > tmp1
cat TSS_hg38.bed | cut -f 6 | sort | uniq -c | awk '($1<2)' | awk '{print $2}' | sort > tmp2
join tmp2 tmp1 | sed s/" "/\\t/g | sort  > tmpA
join legal_genes.txt tmpA | sed s/" "/\\t/g |  cut -f 2,3,4,5,6,7 > tmpB

sort-bed tmpB  > TSS_gene.txt
bedtools window -a TAD.sorted.bed -b TSS_gene.txt -w 0 > test
cut -f 4 test | sort | uniq -c | awk '($1>=2)' | awk  '{print $2}' | sed s/_/\\t/g | awk '{print $1"_"$2"\t"$2}' | sort -nk2 | sed s/\\t/\;/g > candidate_TAD_list


sh for_loop.sh | awk '($1!=$2)' | sort | uniq | awk '{if($3==$4){print $0"\tsame_direction"}else{print $0"\topposite_direction"}}' > Gene_pairs_with_sampe_TAD.txt 

cat Gene_pairs_with_sampe_TAD.txt | awk '{if($1<$2){print $1"\t"$2"\t"$5}else{print $2"\t"$1"\t"$5}}' | sort | uniq > uniq_Gene_pairs_with_sampe_TAD.txt 

cat TSS_gene.txt  | awk '{print $6"\t"($2+$3)/2}' | sort > gene_position.txt

join uniq_Gene_pairs_with_sampe_TAD.txt gene_position.txt | sed s/" "/\\t/g | sort -k2 > sortedLeft
join sortedLeft gene_position.txt -1 2 -2 1 | sed s/" "/\\t/g | sort | awk '{if($4>$5){print $1"\t"$2"\t"$3"\t"$4-$5}else{print $1"\t"$2"\t"$3"\t"$5-$4}}' |  awk '{if($1<$2){print $0}else{print $2"\t"$1"\t"$3"\t"$4}}' > ../final.TAD_gene_pairs.txt

