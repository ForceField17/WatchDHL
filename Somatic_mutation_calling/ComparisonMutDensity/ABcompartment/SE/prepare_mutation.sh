tail -n+2 ../../Final_results/DHL.specific.seg | awk '{print $2"\t"$3"\t"$4"\t"$1}' > tmp 
sort-bed tmp > DHL

tail -n+2 ../../Final_results/FL.specific.seg | awk '{print $2"\t"$3"\t"$4"\t"$1}' > tmp
sort-bed tmp | sed s/P8.a/P8/g  > FL
