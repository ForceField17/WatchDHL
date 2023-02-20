tail -n+2 ../../Final_results/DHL.specific.seg | awk '{if(length($5)>1 || length($6)>1){print $0"\tIndel"}else{print $0"\tSNV"}}' | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$1"\t"$7"\t"$8}' > tmp
sort-bed tmp > theDHL
tail -n+2 ../../Final_results/FL.specific.seg | awk '{if(length($5)>1 || length($6)>1){print $0"\tIndel"}else{print $0"\tSNV"}}' | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$1"\t"$7"\t"$8}' > tmp
sort-bed tmp > theFL

perl anno_signature.pl theFL > FL_AID &
perl anno_signature.pl theDHL > DHL_AID &
wait
