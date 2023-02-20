cat FL_somatic_Mut.specific.txt | awk '($105!="SNV"  )' | cut -f 1,2,3,4,5 > exist_FL_indel.txt
cat DHL_somatic_Mut.specific.txt | awk '($105!="SNV" )' | cut -f 1,2,3,4,5 > exist_DHL_indel.txt

cut -f 1,2,3,4,5 indel_DHL_somatic_Mut.specific.txt > tmp_DHL
cut -f 1,2,3,4,5 indel_DHL_somatic_Mut.specific.txt >> tmp_DHL

cut -f 1,2,3,4,5 indel_FL_somatic_Mut.specific.txt >  tmp_FL
cut -f 1,2,3,4,5 indel_FL_somatic_Mut.specific.txt >> tmp_FL

cat exist_FL_indel.txt tmp_FL | sort | uniq -c | sort | awk '{print $1"\t"$0}' > the_FL
cat exist_DHL_indel.txt tmp_DHL | sort | uniq -c | sort  | awk '{print $1"\t"$0}' > the_DHL

rm tmp_*
rm exist_*.txt
