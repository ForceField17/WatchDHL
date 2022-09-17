tail -n+2 ./step4.FL.csv > temp
tail -n+2 ./step4.DHL.csv >> temp
cat temp  | cut -f 19,20,21,117 | awk '($2!=".")' | sort | uniq -c | sort -nr | awk '{ print $3"\t"$1"\t"$2"\t"$4"\t"$5}' | sort | awk -F"\t" '($2>=1)' > the_tmp.txt 
join -a1 -1 1 the_tmp.txt -2 1 ./OncKB.sorted.list | sed s/" "/\\t/g  | awk -F"\t" '{if($6!=""){print $0}else{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t-"}}' | sed s/"Oncogene,TSG"/ambiguous/g | sort > the_tmp2.txt
join -j1 the_tmp2.txt ./known_driver.list -a1 | sed s/" "/\\t/g | awk -F"\t" '{if($7!=""){print $0}else{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t-"}}' | sort -nrk2
