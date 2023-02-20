#!/bin/bash


../savi.env
set -e
export PATH=$PATH:/scratch/PI/jgwang/dsongad/software/snpEff
export LIBRARY_PATH=/scratch/PI/jgwang/dsongad/software/ncurses/lib/:$LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/PI/jgwang/dsongad/software/ncurses/lib/

OUTPUT=$1
SN=$2
LG=$3
HG=$4

V=/scratch/PI/jgwang/dsongad/software/my_lib/tmp_vcf_lib
S=/scratch/PI/jgwang/dsongad/IGCT/run_savi/SAVI
R=/scratch/PI/jgwang/dsongad/software/my_lib/star_genome_d1_vd1_gtfv22/GRCh38.d1.vd1.fa

for i in {1..22} {X,Y}; do
	mkdir -p ${OUTPUT}/${i}
        ${S}/savi_SONG_allregion.py --bams ${SN},${LG},${HG} \
                --names SuNor,LG,HG \
                --ann hg38_Song --memory 12 \
                --ref ${R} -c 2:1,3:1 \
                --region chr${i} -v --superverbose \
                --noncoding --mapqual 20 --mindepth 4 \
                --outputdir ${OUTPUT}/${i} \
                --annvcf ${V}/Cosmic_2020_hg38/Song_hg38_Cosmic.sorted.vcf.gz,${V}/Cosmic_2020_hg38/Song_hg38_ClinVar.sorted.vcf.gz,${V}/Cosmic_2020_hg38/Song_hg38_PancancerGermline.sorted.vcf.gz,${V}/Cosmic_2020_hg38/Song_hg38_cbio.vcf.gz,${V}/normal_hg38/Song_hg38_All_SNP.sorted.vcf.gz,${V}/normal_hg38/Song_hg38_TOPMED.sorted.vcf.gz,${V}/normal_hg38/IGCTfamB.sorted.vcf.gz,${V}/normal_hg38/IGCTpatB.sorted.vcf.gz,${V}/normal_hg38/Song_hg38_GATKnormal.vcf.gz,${V}/normal_hg38/Song_hg38_HGDP_1KG.sorted.vcf.gz,${V}/normal_hg38/Song_hg38_MuTectPON.vcf.gz,${V}/normal_hg38/Song_hg38_dbSNP_ALFA.sorted.vcf.gz,${V}/normal_hg38/Song_hg38_gnomAD.sorted.vcf.gz,${V}/normal_hg38/Song_hg38_219normals.vcf.gz,${V}/normal_hg38/Song_China_map.sorted.vcf.gz,${V}/normal_hg38/Song_hg38_ExAC.sorted.vcf.gz,${V}/normal_hg38/Song_hg38_NGDC_SNP.vcf.gz,${V}/normal_hg38/Song_hg38_meganormal.vcf.gz 1> ${OUTPUT}/${i}/err.out 2>${OUTPUT}/${i}/err.log &
done
wait

echo SAVI done!

echo Start merging...

O1=${OUTPUT}/merge.report.coding.PDfilter.txt
O2=${OUTPUT}/merge.report.coding.txt


head -n 1 ${OUTPUT}/1/report/report.coding.PDfilter.txt > ${O1}
head -n 1 ${OUTPUT}/1/report/report.coding.txt > ${O2}

for i in {1..22} {X,Y}; do
        F1=${OUTPUT}/${i}/report/report.coding.PDfilter.txt
	F2=${OUTPUT}/${i}/report/report.coding.txt
	
        if [ -e ${F1} ]; then
                tail -n +2 ${F1} >> ${O1}
        else
                echo ${F1} does not exist
        fi

        if [ -e ${F2} ]; then
                tail -n +2 ${F2} >> ${O2}
        else
                echo ${F2} does not exist
        fi

done

