NNN=$1

XXX=$2

mkdir TMP.$NNN

cat $NNN > TMP.$NNN/MYFILE
#head -n1 $NNN | awk '{print $0"\tNormalMQ0"}' > TMP.$NNN/header.txt
cd TMP.$NNN
a=(`wc -l MYFILE`) ; lines=`echo $(($a/40))` ; split -l $lines -d  MYFILE file

for i in file*;do
	cut -f 2,3 $i > bed_$i\.bed
	nohup sambamba-0.8.0-linux-amd64-static -q depth base --nthreads=1 -L ./bed_$i\.bed -z ../../../../../../raw_data/results_B_DNA/merged_Normal.bam -F 'mapping_quality == 0' -o zzz_$i &
done

wait

for i in file*;do
	tail -n+2 zzz_$i | awk -F"\t" '{print $1"_"$2"\t"$3}' > xxx_$i
	cat bed_$i\.bed | awk '{count++} {print $1"_"$2"\t"count}' | sort > tmp1_$i
	cat xxx_$i | sort > tmp2_$i
	join tmp1_$i tmp2_$i -a1 | sort -nk2 | awk '{print $3}' > yyy_$i
	paste $i yyy_$i > NEW_$i
done
wait


cat NEW*  |  sed s/\\r//g > ../$XXX
cd ..
rm -r TMP.$NNN
