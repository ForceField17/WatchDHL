cat  ~/gencode.v40.basic.annotation.gtf | awk '{if($3 == "gene") {if($7 == "-") {$4=$5-50 ; $5+=50} else {$5=$4+50 ; $4-=50} print}}' |sed s/" "/\\t/g |  cut -f 1,4,5,7,12,14 | sed s/\"//g | sed s/\;//g  | awk -F"\t" '($5~/^IG/ || $5~/^protein/ || $5~/^lncRNA/)' > TSS_hg38.bed

cut -f 1,2,3 TSS_hg38.bed > tmp
sort-bed tmp > xxx
bedtools merge -i xxx > Input_peaks.bed
