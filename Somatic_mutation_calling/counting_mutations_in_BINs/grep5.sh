NNN=$1
FILE=$2
cat $FILE | awk -v XXX=$NNN '($120==XXX)' | cut -f 124,125,128 | tr '\t' '\n' |  sort | uniq | grep -v "\." |  tr '\n' ','
