#!/bin/bash

rawName=$1
cleanName=$2

fastp -w 12 -i ../$rawName\.R1.fastq.gz -o ./$cleanName\.R1.fastq.gz \
 -I ../$rawName\.R2.fastq.gz -O ./$cleanName\.R2.fastq.gz \
 -q 15 -u 40 --length_required 45 \
 --detect_adapter_for_pe \
 -h fastp_$cleanName\.html -j fastp_$cleanName\.json -R report_$cleanName\.txt
