cat ~/hg38_clean.basic.txt | awk '($7=="gene" && $6!~/^ENSG[0-9]+/)' | awk '($6~/^MIR[0-9]/  || $5=="protein_coding")' | cut -f 1,2,3,4,5,6 > gene.bed 
echo -en "chr14\t105586437\t106879844\t.\t.\tIGH\t.\n" >> gene.bed
echo -en "chr22\t22026076\t22922913\t.\t.\tIGL\t.\n" >> gene.bed
echo -en "chr2\t88857361\t90235368\t.\t.\tIGK\t.\n" >> gene.bed
echo -en "chr3\t187889202\t188007603\t.\t.\tBCL6-SE\tSE\n" >> gene.bed

cat ../step4.*.csv  | awk '($4~/P[126789]/ && ($14~/^int/ || $14~/UTR/) )' |  cut -f 20 | awk '($1!="" && $1!="-" && $1!=".")' | sort | uniq -c | awk '($1>=2)'  | awk '{print $2}' | tr "," "\n" | sort > b
cat ../ranking_SE_Histone_3kb.txt  | awk '($5>=5 && $8>=2 && $11!="")'  | cut -f 11 | awk '($1!="")' | sort | uniq -c | sort | awk '{print $2}' | tr "," "\n" | sort > c
join b c > a
echo "MIR142" >> a
echo "BCL6-SE" >> a
awk 'NR==FNR {a[$1]++; next} $6 in a' a gene.bed  > b


sort-bed b | awk '{print $1"\t"$2"\t"$3"\t"$6}' | awk '{print $0"\t"$3-$2}' | grep -v chrM  | awk '($2>0)' > sorted.bed

tail -n+2 ../step4.FL.csv  | awk '($4~/P[126789]/ && ($14~/^int/ || $14~/UTR/) )' | cut -f 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26,35,38,51,53,54,55,75,84,93  > FL.bed
tail -n+2 ../step4.DHL.csv | awk '($4~/P[126789]/ && ($14~/^int/ || $14~/UTR/) )' | cut -f 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26,35,38,51,53,54,55,75,84,93 > DHL.bed


bedtools closest -a FL.bed -b sorted.bed -t first -d | awk -F"\t" '($40<=50000 && $40!=-1)' > annotated.FL
bedtools closest -a DHL.bed -b sorted.bed -t first -d  | awk -F"\t" '($40<=50000 && $40!=-1)' > annotated.DHL

cut -f 38,39 annotated.DHL annotated.FL  | sort | uniq -c | awk '($1>=5)' | awk '($1/$3*1000000>=50 || $1>=10)' | awk '{print $2"\t"$1}' | sort > AAA
cut -f 4,38 annotated.DHL annotated.FL | sort | uniq | awk '{print $2"\t"$1}' | sort > BBB
join AAA BBB | sed s/" "/\\t/g | cut -f 1,2 | sort | uniq -c |  awk '{print $2"\t"$3"\t"$1}' | sort > CCC

cut -f 38,39 annotated.FL  | sort | uniq -c | awk '{print $2"\t"$1}' | sort > AAA
cut -f 4,38  annotated.FL  | sort | uniq    | awk '{print $2"\t"$1}' | sort > BBB
join AAA BBB | sed s/" "/\\t/g | cut -f 1,2 | sort | uniq -c |  awk '{print $2"\t"$3"\t"$1}' | sort > DDD


cut -f 38,39 annotated.DHL  | sort | uniq -c | awk '{print $2"\t"$1}' | sort > AAA
cut -f 4,38  annotated.DHL  | sort | uniq    | awk '{print $2"\t"$1}' | sort > BBB
join AAA BBB | sed s/" "/\\t/g | cut -f 1,2 | sort | uniq -c |  awk '{print $2"\t"$3"\t"$1}' | sort > EEE

join CCC DDD | sed s/" "/\\t/g | sort > FFF
echo -en "gene\tTotalNum\tTotalPat\tFLNum\tFLPat\tDHLNum\tDHLPat\trank\n" > number_summary.txt
join FFF EEE | sed s/" "/\\t/g | sort -nk2 -nk3 | awk '($7>=2 && $6>=5)' | sort -nrk6 | awk '{count++} {print $0"\t"count}' >> number_summary.txt





head -n1 ../step4.FL.csv  | cut -f 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26,35,38,51,53,54,55,75,84,93 | awk '{print $0"\ta\tb\tc\tGeneID\tsize\tdistance"}' > annotated.frequent.txt
cat annotated.DHL annotated.FL | sort | uniq >> annotated.frequent.txt

chmod 755 number_summary.txt
chmod 755 annotated.frequent.txt
