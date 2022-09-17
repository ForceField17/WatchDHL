
cat ../step4.DHL.csv | awk -F"\t" '($15!="LOW" && $15!="MODIFIER")' > step4.DHL.csv
cat ../step4.FL.csv | awk -F"\t" '($15!="LOW" && $15!="MODIFIER")' > step4.FL.csv
perl countES.pl ranks > ranking_gene_frequency_table.txt
wait
