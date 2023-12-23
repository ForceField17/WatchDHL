head -n1 results_table.test.txt | awk '{print $0"\tMutation"}' > gene_pairs_results_table.txt

tail -n+2 results_table.test.txt | cut -f 1,2,3,4,5,6,7 | sort > AAA
cat hyper_genes.txt | awk '{print $1"\thypermutated"}' | sort > CCC
join -a 1 AAA CCC | tr " " "\t"  | awk '{if($8==""){print $0"\tWT"}else{print $0}}' | awk '{print $2"\t"$1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' | sort  > BBB
join -a 1 BBB CCC |  tr " " "\t"  | awk '{if($9==""){print $0"\tWT"}else{print $0}}' | awk '{if($8=="WT" && $9=="WT"){print $2"\t"$1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\tWT"}else{print $2"\t"$1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\tHypermutated"}}' | sort -nk5 >> gene_pairs_results_table.txt

rm AAA BBB CCC

 
