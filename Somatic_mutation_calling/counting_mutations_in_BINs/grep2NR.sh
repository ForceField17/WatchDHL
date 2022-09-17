NNN=$1
FILE=$2
cat $FILE | awk -v XXX=$NNN '($120==XXX)' | wc -l
