TMP="TMP"

cat ../../annotation/kataegis/hypermutated_loci.bed | cut -f 1,2,3 | awk '{print $0"\tB_HyperM"}' > Input_peaks.bed

theL=`cat Input_peaks.bed | wc -l`
theN=`cat Input_peaks.bed | awk '{print $3-$2}' | awk '{x++;qq+=$1}END{print int(qq/x)}'`
cat ../../annotation/generate_BINs/contig.bed | grep -v chrM | cut -f 1,3 > hg38.genome
bedtools random -l $theL -n $theN -g hg38.genome | awk '{print $1"\t"$2"\t"$3"\tA_Random"}' >> Input_peaks.bed

cat ../../annotation/generate_BINs/contig.bed | grep -v chrM > xxx
bedtools subtract -a xxx -b ../../annotation/kataegis/hypermutated_loci.bed | awk '{print $1"\t"$2"\t"$3"\tA_Background"}' >>  Input_peaks.bed

sort-bed Input_peaks.bed > AAA
mv AAA Input_peaks.bed

cat Input_peaks.bed | awk '($4=="B_HyperM"                    )' | awk '{print $3-$2}' | awk '{count+=$1} END {print count}'  | awk '{print "B_HyperM"                   , $1}' > $TMP\3
cat Input_peaks.bed | awk '($4=="A_Random"                     )' | awk '{print $3-$2}' | awk '{count+=$1} END {print count}'  | awk '{print "A_Random"                    , $1}' >> $TMP\3
cat Input_peaks.bed | awk '($4=="A_Background"                     )' | awk '{print $3-$2}' | awk '{count+=$1} END {print count}'  | awk '{print "A_Background"                    , $1}' >> $TMP\3

cat $TMP\3 | sed s/" "/\\t/g | sort > $TMP\4

echo -en "xxx\ttags\tcase\ttumor\tcompartment\tmutation\ttotal_mut\tsize\n"

bedtools window -a Input_peaks.bed -b nonAID/FL -w 0 | cut -f 4,8 | sort | uniq -c | awk '{print $3"+"$2"\t"$1}' | sort > $TMP\1
join selectedFL $TMP\1 -a1 | sed s/" "/\\t/g | awk '{if($2==""){print $1"\t0"}else{print $0}}' | sed s/\+/\\t/g > $TMP\6
cut -f 4 FL | sort | uniq -c | awk '{print $2"\t"$1}' | sort > $TMP\2
join $TMP\6 $TMP\2 | sed s/" "/\\t/g | awk '{print $2"\t"$0}' | sort > $TMP\5
join $TMP\5 $TMP\4 | sed s/" "/\\t/g | cut -f 2,3,4,5,6,7 | sed s/.FL/\\tFL/g | awk '{print "X_"$3"\t"$1"."$2"\t"$0}'


bedtools window -a Input_peaks.bed -b nonAID/DHL -w 0 | cut -f 4,8 | sort | uniq -c | awk '{print $3"+"$2"\t"$1}' | sort > $TMP\1
join selectedDHL $TMP\1 -a1 | sed s/" "/\\t/g | awk '{if($2==""){print $1"\t0"}else{print $0}}' | sed s/\+/\\t/g > $TMP\6
cut -f 4 DHL | sort | uniq -c | awk '{print $2"\t"$1}' | sort > $TMP\2
join $TMP\6 $TMP\2 | sed s/" "/\\t/g | awk '{print $2"\t"$0}' | sort > $TMP\5
join $TMP\5 $TMP\4 | sed s/" "/\\t/g | cut -f 2,3,4,5,6,7 | sed s/.DHL/\\tDHL/g | awk '{print "Y_"$3"\t"$1"."$2"\t"$0}'

