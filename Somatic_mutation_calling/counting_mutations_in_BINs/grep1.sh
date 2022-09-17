tail -n+2 ./step4.FL.csv > temp
tail -n+2 ./step4.DHL.csv >> temp
cat temp  | cut -f 117,118,119,120 | awk '($4!=".")' | sort | uniq -c | sort -nr | awk '{ print $5"\t"$1"\t"$2"\t"$3"\t"$4}' | sort -nrk2 | awk -F"\t" '($2>=1)'

