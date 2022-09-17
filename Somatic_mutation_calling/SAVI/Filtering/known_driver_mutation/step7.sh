 
	cat	../annotated.*.txt  | awk '{print $1"_"$2"_"$3"_"$4"_"$5"\t"$0}' | sort | uniq > BBB
	join AAA BBB -a1 | sed s/" "/\\t/g > annotated.txt &
wait

        cat     ../indel_PN_annotated.*.txt  | awk '{print $1"_"$2"_"$3"_"$4"_"$5"\t"$0}' | sort | uniq > BBB
	join CCC BBB -a1 | sed s/" "/\\t/g > indel_annotated.txt &

wait
