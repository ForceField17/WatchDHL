
cat ../OCI.Ly1.H3K4me3.normal.ENCFF670JNI.bed | awk '{print $4"\t"$1":"$2"-"$3"\t"$3-$2"\t"$5"\t"$7"\t"$8}' | sed s/Peak_//g | sort > tmp_H3K4me3

echo -en "Rank\tPos\tsize\tscore\tEnrichment\tpvalue\tMutNumber\tTumor\n" >  Table_H3K4me3.PeakMut.txt
bedtools window -b DHL -a ../OCI.Ly1.H3K4me3.normal.ENCFF670JNI.bed -w 0 | cut -f 1,2,3,4,5,7,8 | sort | uniq -c | awk '{print $5"\t"$1}' | sed s/Peak_//g | sort > tmp_DHL 
join tmp_H3K4me3 tmp_DHL -a 1 | sed s/" "/\\t/g | awk '{if($7==""){print $0"\t0\tDHL"}else{print $0"\tDHL"}}' >> Table_H3K4me3.PeakMut.txt

bedtools window -b FL -a ../OCI.Ly1.H3K4me3.normal.ENCFF670JNI.bed -w 0 | cut -f 1,2,3,4,5,7,8 | sort | uniq -c | awk '{print $5"\t"$1}' | sed s/Peak_//g | sort > tmp_FL
join tmp_H3K4me3 tmp_FL -a 1 | sed s/" "/\\t/g | awk '{if($7==""){print $0"\t0\tFL"}else{print $0"\tFL"}}' >> Table_H3K4me3.PeakMut.txt


cat ../OCI.Ly1.H3K27ac.normal.ENCFF909YRO.bed  | awk '{print $4"\t"$1":"$2"-"$3"\t"$3-$2"\t"$5"\t"$7"\t"$8}' | sed s/Peak_//g | sort > tmp_H3K27ac

echo -en "Rank\tPos\tsize\tscore\tEnrichment\tpvalue\tMutNumber\tTumor\n" >  Table_H3K27ac.PeakMut.txt
bedtools window -b DHL -a ../OCI.Ly1.H3K27ac.normal.ENCFF909YRO.bed -w 0 | cut -f 1,2,3,4,5,7,8 | sort | uniq -c | awk '{print $5"\t"$1}' | sed s/Peak_//g | sort > tmp_DHL
join tmp_H3K27ac tmp_DHL -a 1 | sed s/" "/\\t/g | awk '{if($7==""){print $0"\t0\tDHL"}else{print $0"\tDHL"}}' >> Table_H3K27ac.PeakMut.txt

bedtools window -b FL -a ../OCI.Ly1.H3K27ac.normal.ENCFF909YRO.bed -w 0 | cut -f 1,2,3,4,5,7,8 | sort | uniq -c | awk '{print $5"\t"$1}' | sed s/Peak_//g | sort > tmp_FL
join tmp_H3K27ac tmp_FL -a 1 | sed s/" "/\\t/g | awk '{if($7==""){print $0"\t0\tFL"}else{print $0"\tFL"}}' >> Table_H3K27ac.PeakMut.txt



cat ../OCI.Ly1.H3K4me1.normal.ENCFF348SPV.bed   | awk '{print $4"\t"$1":"$2"-"$3"\t"$3-$2"\t"$5"\t"$7"\t"$8}' | sed s/Peak_//g | sort > tmp_H3K4me1

echo -en "Rank\tPos\tsize\tscore\tEnrichment\tpvalue\tMutNumber\tTumor\n" >  Table_H3K4me1.PeakMut.txt
bedtools window -b DHL -a ../OCI.Ly1.H3K4me1.normal.ENCFF348SPV.bed  -w 0 | cut -f 1,2,3,4,5,7,8 | sort | uniq -c | awk '{print $5"\t"$1}' | sed s/Peak_//g | sort > tmp_DHL
join tmp_H3K4me1 tmp_DHL -a 1 | sed s/" "/\\t/g | awk '{if($7==""){print $0"\t0\tDHL"}else{print $0"\tDHL"}}' >> Table_H3K4me1.PeakMut.txt

bedtools window -b FL -a ../OCI.Ly1.H3K4me1.normal.ENCFF348SPV.bed  -w 0 | cut -f 1,2,3,4,5,7,8 | sort | uniq -c | awk '{print $5"\t"$1}' | sed s/Peak_//g | sort > tmp_FL
join tmp_H3K4me1 tmp_FL -a 1 | sed s/" "/\\t/g | awk '{if($7==""){print $0"\t0\tFL"}else{print $0"\tFL"}}' >> Table_H3K4me1.PeakMut.txt


