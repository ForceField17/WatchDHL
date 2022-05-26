#!/usr/bin/bash
FUSION=$1
FIGURE=$2
BAM=$3

draw_fusions.R --annotation=./Gencode.v22.SVgenes.gtf --fusions=$FUSION --output=$FIGURE --cytobands=/home/dsongad/software/arriba_v2.0.0/database/cytobands_hg38_GRCh38_v2.0.0.tsv --proteinDomains=~/software/arriba_v2.0.0/database/protein_domains_hg38_GRCh38_v2.0.0.gff3

