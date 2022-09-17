NNN=$1
FILE=$2

cat $FILE | awk -v XXX=$NNN '($20==XXX)' | cut -f 4,15 |  sort | cut -f 2 | grep -v "\." |  tr '\n' ','
