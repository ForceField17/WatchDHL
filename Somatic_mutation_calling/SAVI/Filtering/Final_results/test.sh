CCC=$1
 AAA=$2

BBB=$3
 samtools view -b /scratch/PI/jgwang/dsongad/DLBCL/raw_data/results_T_DNA/P7_LG.sorted.MD.bam $CCC\:$AAA\-$AAA | samtools calmd -e - /scratch/PI/jgwang/dsongad/software/my_lib/star_genome_d1_vd1_gtfv22/GRCh38.d1.vd1.fa | grep -v "^@" | awk -F"\t" '($5>=30)' > tmp


 cat tmp | awk -v pos=$AAA  'BEGIN {OFS = FS = "\t" } ; {n=split($10,a,"") ; if(a[(pos-$4)+1] == "=" ) print pos-$4+1 "\t" length($10)-pos+$4-1 "\t" $12 "\t" $10}' | awk '($1>=37 && $2>=37 && ($1+$2)>=113 && $3!~/^XA/)' | sort | uniq 

cat tmp | awk -v pos=$AAA -v theB=$BBB 'BEGIN {OFS = FS = "\t" } ; {n=split($10,a,"") ; if(a[(pos-$4)+1] == theB ) print pos-$4+1 "\t" length($10)-pos+$4-1 "\t" $12 "\t" $10}' | awk '($1>=37 && $2>=37 && ($1+$2)>=113 && $3!~/^XA/)' | sort | uniq | cut -f 4 | sed s/=//g | awk '{print length($1)}'

#cat tmp | awk -v pos=$AAA -v theB=$BBB 'BEGIN {OFS = FS = "\t" } ; {n=split($10,a,"") ; if(a[(pos-$4)+1] == theB ) print pos-$4+1 "\t" length($10)-pos+$4-1 "\t" $12 "\t" $10}' #| awk '($1>=37 && $2>=37 && ($1+$2)>=113 && $3!~/^XA/)' | sort | uniq
