TMP="TMP"

cat ../../DistanceMutDensity/Slide_Both_H3K4me3_and_H3K27ac/Input_peaks.bed | cut -f 1,2,3 | awk '{print $0"\tH_Intersection_4me3_and_27ac"}' > Input_peaks.bed
cat ../../DistanceMutDensity/Slide_H3K27ac/Input_peaks.bed | cut -f 1,2,3 | awk '{print $0"\tB_H3K27ac"}' >> Input_peaks.bed
cat ../../DistanceMutDensity/Slide_H3K4me3/Input_peaks.bed | cut -f 1,2,3 | awk '{print $0"\tC_H3K4me3"}' >> Input_peaks.bed
cat ../../DistanceMutDensity/Slide_H3K4me1/Input_peaks.bed | cut -f 1,2,3 | awk '{print $0"\tD_H3K4me1"}' >> Input_peaks.bed
cat ../../DistanceMutDensity/Slide_Only_H3K27ac/Input_peaks.bed | cut -f 1,2,3 | awk '{print $0"\tF_Only_27ac"}' >> Input_peaks.bed
cat ../../DistanceMutDensity/Slide_Only_H3K4me3/Input_peaks.bed | cut -f 1,2,3 | awk '{print $0"\tG_Only_4me3"}' >> Input_peaks.bed



theL=`cat Input_peaks.bed | awk '($4~/H3K/)' | wc -l | awk '{print int($1/3)}'`
theN=`cat Input_peaks.bed | awk '($4~/H3K/)'| awk '{print $3-$2}' | awk '{x++;qq+=$1}END{print int(qq/x)}'`
cat ../../annotation/generate_BINs/contig.bed | grep -v chrM | cut -f 1,3 > hg38.genome
bedtools random -l $theL -n $theN -g hg38.genome | awk '{print $1"\t"$2"\t"$3"\tA_Random"}' >> Input_peaks.bed


theL=`cat Input_peaks.bed | awk '($4~/Only/ || $4~/Inter/)' | wc -l | awk '{print int($1/3)}'`
theN=`cat Input_peaks.bed | awk '($4~/Only/ || $4~/Inter/)' | awk '{print $3-$2}' | awk '{x++;qq+=$1}END{print int(qq/x)}'`
cat ../../annotation/generate_BINs/contig.bed | grep -v chrM | cut -f 1,3 > hg38.genome
bedtools random -l $theL -n $theN -g hg38.genome | awk '{print $1"\t"$2"\t"$3"\tE_Random"}' >> Input_peaks.bed


sort-bed Input_peaks.bed > AAA
mv AAA Input_peaks.bed


cat Input_peaks.bed | awk '($4=="B_H3K27ac"                    )' | awk '{print $3-$2}' | awk '{count+=$1} END {print count}'  | awk '{print "B_H3K27ac"                   , $1}' > $TMP\3
cat Input_peaks.bed | awk '($4=="C_H3K4me3"                    )' | awk '{print $3-$2}' | awk '{count+=$1} END {print count}'  | awk '{print "C_H3K4me3"                   , $1}' >> $TMP\3
cat Input_peaks.bed | awk '($4=="D_H3K4me1"                    )' | awk '{print $3-$2}' | awk '{count+=$1} END {print count}'  | awk '{print "D_H3K4me1"                   , $1}' >> $TMP\3
cat Input_peaks.bed | awk '($4=="F_Only_27ac"                  )' | awk '{print $3-$2}' | awk '{count+=$1} END {print count}'  | awk '{print "F_Only_27ac"                 , $1}' >> $TMP\3
cat Input_peaks.bed | awk '($4=="G_Only_4me3"                  )' | awk '{print $3-$2}' | awk '{count+=$1} END {print count}'  | awk '{print "G_Only_4me3"                 , $1}' >> $TMP\3
cat Input_peaks.bed | awk '($4=="H_Intersection_4me3_and_27ac" )' | awk '{print $3-$2}' | awk '{count+=$1} END {print count}'  | awk '{print "H_Intersection_4me3_and_27ac", $1}' >> $TMP\3
cat Input_peaks.bed | awk '($4=="A_Random"                     )' | awk '{print $3-$2}' | awk '{count+=$1} END {print count}'  | awk '{print "A_Random"                    , $1}' >> $TMP\3
cat Input_peaks.bed | awk '($4=="E_Random"                     )' | awk '{print $3-$2}' | awk '{count+=$1} END {print count}'  | awk '{print "E_Random"                    , $1}' >> $TMP\3

cat $TMP\3 | sed s/" "/\\t/g | sort > $TMP\4

echo -en "case\ttumor\tcompartment\tmutation\ttotal_mut\tsize\n"

bedtools window -a Input_peaks.bed -b FL -w 0 | cut -f 4,8 | sort | uniq -c | awk '{print $3"+"$2"\t"$1}' | sort > $TMP\1
join selectedFL $TMP\1 -a1 | sed s/" "/\\t/g | awk '{if($2==""){print $1"\t0"}else{print $0}}' | sed s/\+/\\t/g > $TMP\6
cut -f 4 FL | sort | uniq -c | awk '{print $2"\t"$1}' | sort > $TMP\2
join $TMP\6 $TMP\2 | sed s/" "/\\t/g | awk '{print $2"\t"$0}' | sort > $TMP\5
join $TMP\5 $TMP\4 | sed s/" "/\\t/g | cut -f 2,3,4,5,6,7 | sed s/.FL/\\tFL/g


bedtools window -a Input_peaks.bed -b DHL -w 0 | cut -f 4,8 | sort | uniq -c | awk '{print $3"+"$2"\t"$1}' | sort > $TMP\1
join selectedDHL $TMP\1 -a1 | sed s/" "/\\t/g | awk '{if($2==""){print $1"\t0"}else{print $0}}' | sed s/\+/\\t/g > $TMP\6
cut -f 4 DHL | sort | uniq -c | awk '{print $2"\t"$1}' | sort > $TMP\2
join $TMP\6 $TMP\2 | sed s/" "/\\t/g | awk '{print $2"\t"$0}' | sort > $TMP\5
join $TMP\5 $TMP\4 | sed s/" "/\\t/g | cut -f 2,3,4,5,6,7 | sed s/.DHL/\\tDHL/g

