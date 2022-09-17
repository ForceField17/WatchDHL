NNN=$1
FILE=$2

cat $FILE | awk -v XXX=$NNN '($120==XXX)' | cut -f 122,129 | tr '\t' '\n' |  sort | uniq | grep -v "\." |  tr '\n' ','
