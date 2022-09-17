TUMOR=$1
TMP=$2
TOTAL=`cat $TUMOR | cut -f 4 | sort | uniq -c | awk '{print $2"\t"$1}' | cut -f 2 | tr "\n" ","`
cp Input_peaks.bed $TMP\_old

SIZE=`cat Input_peaks.bed | awk '{print $3-$2}' | awk '{count+=$1} END {print count}'`
bedtools window -a Input_peaks.bed -b $TUMOR -w 0 | cut -f 8 | sort | uniq -c | awk '{print $2"\t"$1}' | sort | sed s/.DHL//g | sed s/.FL//g > $TMP\4
NUM=`join selected $TMP\4 -a1 | awk '{if($2==""){print 0}else{print $2}}' | tr "\n" ","`
echo -en $TUMOR"\t0\t0\t"$NUM""$SIZE"\t"$TOTAL"\n"
j=0
#for i in {100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,4000,5000,6000,7000,8000,9000,10000,12000,14000,16000,18000,20000,25000,30000,40000,50000,60000,80000,100000,150000,200000,300000,400000,800000};do
#for i in {1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000,16000,17000,18000,19000,20000,21000,22000,23000,24000,25000,26000,27000,28000,29000,30000,31000,32000,33000,34000,35000,36000,37000,38000,39000,40000,41000,42000,43000,44000,45000,46000,47000,48000,49000,50000};do
for n in `cat list1`;do
    i=`echo $n | cut -d "," -f 1`
    j=`echo $n | cut -d "," -f 2`
cat Input_peaks.bed | awk -F"\t" -v DIS=$i '{if($4~/^-/){print $1"\t"$2+DIS"\t"$3+DIS}else{print $1"\t"$2-DIS"\t"$3-DIS}}' |  awk -F"\t" '{if($2<0){print $1"\t0\t"$3}else{print $0}}' |  awk -F"\t" '{if($3<0){print $1"\t"$2"\t1"}else{print $0}}' > $TMP\0
	sort-bed $TMP\0 >$TMP\1	
	cat Input_peaks.bed | awk -F"\t" -v DIS=$i '{if($4~/^-/){print $1"\t"$2"\t"$2+DIS}else{print $1"\t"$3-DIS"\t"$3}}' |  awk -F"\t" '{if($2<0){print $1"\t0\t"$3}else{print $0}}' |  awk -F"\t" '{if($3<0){print $1"\t"$2"\t1"}else{print $0}}' > $TMP\_old
	bedtools subtract -a $TMP\1 -b $TMP\_old > $TMP\2_1
	sort-bed $TMP\2_1 > $TMP\2
	bedtools merge -i $TMP\1 > $TMP\_old
	bedtools merge -i $TMP\2 > $TMP\3
	SIZE=`cat $TMP\3 | awk '{print $3-$2}' | awk '{count+=$1} END {print count}'`
	bedtools window -a $TMP\3 -b $TUMOR -w 0 | cut -f 7 | sort | uniq -c | awk '{print $2"\t"$1}' | sort | sed s/.DHL//g | sed s/.FL//g > $TMP\4
        NUM=`join selected $TMP\4 -a1 | awk '{if($2==""){print 0}else{print $2}}' | tr "\n" ","`

	
	echo -en $TUMOR"\t-"$i"\t-"$i"\t"$NUM""$SIZE"\t"$TOTAL"\n"
	j=$i
done	
