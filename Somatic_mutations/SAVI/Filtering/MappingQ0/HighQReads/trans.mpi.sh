NNN=$1
XXX=$2


sh mutType.sh $NNN\.csv > tmp1.$NNN.tmp
sh parallel.sh tmp1.$NNN.tmp $XXX step6.$NNN\.txt



#cp step6.$NNN\.txt ../step6_rmSecondaryAnno/
#cp step6.$NNN\.cSNP.txt ../step6_rmSecondaryAnno/

rm tmp1.$NNN\.tmp
