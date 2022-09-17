NNN=$1
FILE=$2

cat $FILE | awk -v XXX=$NNN '($120==XXX)' | cut -f 20 | grep -v "-" | sort | uniq | tr '\n' ','
