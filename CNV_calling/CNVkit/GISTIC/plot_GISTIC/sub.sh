perl gistic_2_dataframe_for_plotting.pl ../Gistic_DHL8/scores.gistic ../Gistic_DHL8/del_genes.conf_99.txt ../Gistic_DHL8/amp_genes.conf_99.txt > DHL_gistic_table.tsv 
perl gistic_2_dataframe_for_plotting.pl ../Gistic_FL9/scores.gistic ../Gistic_FL9/del_genes.conf_99.txt ../Gistic_FL9/amp_genes.conf_99.txt > FL_gistic_table.tsv 

Rscript plot_gistic_table.r

