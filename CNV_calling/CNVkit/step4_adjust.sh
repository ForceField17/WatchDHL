tail -n+2 CNV_HG.seg >tmpA 
paste CNV_HG.seg tmpA | awk -F"\t" '{if($1==$7 && $2==$8 && ($9-$4)>1 && ($9-$4)<13000 && NR!=1){print $1"\t"$2"\t"$3"\t"$9-1"\t"$5"\t"$6}else{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' | sed s/.cbs.call//g > GISTIC/DHL_8samples.seg

tail -n+2 CNV_LG.seg >tmpA 
paste CNV_LG.seg tmpA | awk -F"\t" '{if($1==$7 && $2==$8 && ($9-$4)>1 && ($9-$4)<13000 && NR!=1){print $1"\t"$2"\t"$3"\t"$9-1"\t"$5"\t"$6}else{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' | sed s/.cbs.call//g > GISTIC/FL_9samples.seg

