cat  ~/gencode.v40.basic.annotation.gtf | awk '{if($3 == "gene") {if($7 == "-") {$4=$5-250 ; $5+=250} else {$5=$4+250 ; $4-=250} print}}' |sed s/" "/\\t/g |  cut -f 1,4,5,7,12,14 | sed s/\"//g | sed s/\;//g  | awk -F"\t" '($5~/^IG/ || $5~/^protein/ || $5~/^lncRNA/)' > TSS_hg38.bed

cut -f 1,2,3,4 TSS_hg38.bed | grep -v chrM > tmp
sort-bed tmp > Input_peaks.bed
