NNN=$1
XXX=$2
ZZZ=$3

mkdir TMP.$NNN

cat $NNN > TMP.$NNN/MYFILE
head -n1 header.txt | awk '{print $0"\tmutType\tRef\tAlt\tseqL\tseqR\tUniqPatern_FL\tHighQread_FL\tUniqPatern_FL2\tHighQread_FL2\tUniqPatern_DHL\tHighQread_DHL"}' > TMP.$NNN/header.txt
cd TMP.$NNN
a=(`wc -l MYFILE`) ; lines=`echo $(($a/40))` ; split -l $lines -d  MYFILE file

for i in file*;do
	perl ../checkReads.mpi.pl $i $XXX > NEW$i &
done
wait

cat header.txt NEW* > ../$ZZZ
cd ..
rm -r TMP.$NNN
