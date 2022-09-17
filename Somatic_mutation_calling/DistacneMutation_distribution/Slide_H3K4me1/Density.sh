TUMOR=$1
TMP=$2
TOTAL=`cat $TUMOR | cut -f 4 | sort | uniq -c | awk '{print $2"\t"$1}' | cut -f 2 | tr "\n" ","`
cp Input_peaks.bed $TMP\_old

SIZE=`cat Input_peaks.bed | awk '{print $3-$2}' | awk '{count+=$1} END {print count}'`
bedtools window -a Input_peaks.bed -b $TUMOR -w 0 | cut -f 8 | sort | uniq -c | awk '{print $2"\t"$1}' | sort | sed s/.DHL//g | sed s/.FL//g > $TMP\4
NUM=`join selected $TMP\4 -a1 | awk '{if($2==""){print 0}else{print $2}}' | tr "\n" ","`
echo -en $TUMOR"\t0\t0\t"$NUM""$SIZE"\t"$TOTAL"\n"

for n in `cat list1`;do
    i=`echo $n | cut -d "," -f 1`
    j=`echo $n | cut -d "," -f 2`
	cat Input_peaks.bed | awk -F"\t" -v DIS=$i '{print $1"\t"$2+DIS"\t"$3+DIS}' |  awk -F"\t" '($2>=0 && $3>$2)' > $TMP\0
	sort-bed $TMP\0 >$TMP\1	
	cat Input_peaks.bed | awk -F"\t" -v DIS=$i '{print $1"\t"$2"\t"$2+DIS}' |  awk -F"\t" '($2>=0 && $3>$2)' > $TMP\_old
	bedtools subtract -a $TMP\1 -b $TMP\_old > $TMP\2a
	sort-bed $TMP\2a > $TMP\2
	bedtools merge -i $TMP\2 > $TMP\3
	SIZE=`cat $TMP\3 | awk '{print $3-$2}' | awk '{count+=$1} END {print count}'`
	bedtools window -a $TMP\3 -b $TUMOR -w 0 | cut -f 7 | sort | uniq -c | awk '{print $2"\t"$1}' | sort | sed s/.DHL//g | sed s/.FL//g > $TMP\4
	NUM=`join selected $TMP\4 -a1 | awk '{if($2==""){print 0}else{print $2}}' | tr "\n" ","`

	echo -en $TUMOR"\t"$i"\t"$i"\t"$NUM""$SIZE"\t"$TOTAL"\n"
	j=$i
done	
