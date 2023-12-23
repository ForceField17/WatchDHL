for n in `cat candidate_TAD_list`;do
    i=`echo $n | cut -d ";" -f 1`
    cat test | awk -v TAD=$i '($4==TAD)' | awk '{print $10";"$8}' > temp1_$i
    for m in `cat temp1_$i`;do
    	for M in `cat temp1_$i`;do
    		j=`echo $m | cut -d ";" -f 1`
    		k=`echo $m | cut -d ";" -f 2`
    		J=`echo $M | cut -d ";" -f 1`
    		K=`echo $M | cut -d ";" -f 2`
    		echo -en $j"\t"$J"\t"$k"\t"$K"\n"
    	done
    done
    rm temp1_$i
done
