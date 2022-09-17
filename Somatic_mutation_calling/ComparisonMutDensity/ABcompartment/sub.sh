TMP="TMP"

SIZE_A=`cat Input_peaks.bed | awk '($4=="A")' | awk '{print $3-$2}' | awk '{count+=$1} END {print count}'`
SIZE_B=`cat Input_peaks.bed | awk '($4=="B")' | awk '{print $3-$2}' | awk '{count+=$1} END {print count}'`



echo -en "case\ttumor\tcompartment\tmutation\ttotal_mut\tsize\n"

bedtools window -a Input_peaks.bed -b FL -w 0 | cut -f 4,8 | sort | uniq -c | awk '{print $3"\t"$2"\t"$1}' | sort > $TMP\1
cut -f 4 FL | sort | uniq -c | awk '{print $2"\t"$1}' | sort > $TMP\2
join $TMP\1 $TMP\2 | sed s/" "/\\t/g | sort | awk -F"\t" -v A=$SIZE_A -v B=$SIZE_B '{if($2=="A"){print $0"\t"A}else{print $0"\t"B}}' | sed s/.FL/\\tFL/g


bedtools window -a Input_peaks.bed -b DHL -w 0 | cut -f 4,8 | sort | uniq -c | awk '{print $3"\t"$2"\t"$1}' | sort > $TMP\1
cut -f 4 DHL | sort | uniq -c | awk '{print $2"\t"$1}' | sort > $TMP\2
join $TMP\1 $TMP\2 | sed s/" "/\\t/g | sort | awk -F"\t" -v A=$SIZE_A -v B=$SIZE_B '{if($2=="A"){print $0"\t"A}else{print $0"\t"B}}' | sed s/.DHL/\\tDHL/g


