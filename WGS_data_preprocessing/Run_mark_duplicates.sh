#!/usr/bin/bash

set -euxo pipefail
export LD_LIBRARY_PATH=/home/qmu/tools/lib:/opt/intel/advisor_xe_2016.1.40.463413/lib64:/home/qmu/tools/xz-5.2.3/tmp/lib 
picard=/home/share/jgwang/softwares/picard/build/libs/picard.jar


NAME=$1

mkdir tmp_$NAME

java -Djava.io.tmpdir=tmp_$NAME \
	-Xmx100g -XX:+UseParallelGC -XX:ParallelGCThreads=4 \
	 -jar ${picard} MarkDuplicates \
	INPUT=${NAME}.sorted.bam \
	OUTPUT=${NAME}.sorted.MD.bam \
	METRICS_FILE=${NAME}.metrics.txt \
	ASSUME_SORTED=true \
	REMOVE_DUPLICATES=true \
	VALIDATION_STRINGENCY=LENIENT

samtools index ${NAME}.sorted.MD.bam

#rm ${NAME}.sorted.bam
#rm ${NAME}.sorted.bam.bai
#rm ${NAME}.sorted.MD.bam.txt
