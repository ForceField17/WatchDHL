NNN=$1
FILE=$2
cat $FILE | awk -v XXX=$NNN '($20==XXX)' | cut -f 4 |  sort | uniq | grep -v "\." |  tr '\n' ','
