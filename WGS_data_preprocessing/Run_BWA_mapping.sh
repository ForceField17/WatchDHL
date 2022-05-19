#!/usr/bin/bash
A='WangLab'
Sample_ID=$1
FAQ1=$2
FAQ2=$3
result=$4


echo [Directory] `pwd`
echo [Machine] `uname -n`
echo [Start] `date`
echo [args] $*
time1=$( date "+%s" )

export PATH=$PATH:/home/share/jgwang/bin/samtools-1.2
export PATH=$PATH:/home/share/jgwang/bin/bwa


bwa mem -t 20 -T 0 -R "@RG\tCN:$A\tID:$Sample_ID" ~/software/my_lib/TCGA_GRCh38_bwa_ref/GRCh38.d1.vd1.fa $FAQ1 $FAQ2 2> tmp.$result | samtools view -Sbh1 -@ 20 - -o $result\.bam

#### mapping duplicates
###counting the num of bam reads

samtools sort -@ 20 -m 2G -O bam -o $result.sorted.bam -T the_temp_$result $result\.bam
samtools index $result\.sorted.bam 
rm $result\.bam

time2=$( date "+%s" )
echo [deltat] $(( $time2 - $time1 ))
echo [End]
