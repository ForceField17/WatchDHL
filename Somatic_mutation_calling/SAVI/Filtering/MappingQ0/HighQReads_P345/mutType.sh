NNN=$1

#echo -en "chr\tstart\tend\tCaseID\tref\talt\trefDepth_B\taltDepth_B\trefDepth_T\taltDepth_T\tAFreq\tGeneName\trefDepth_R\taltDepth_R\tAAchange\tSnpID\tsnpAF\tCosmicID\tnCosmic\tEffect\tImpact\tFunctionClass\tAAlength\tBioType\tGeneCoding\tTranscriptID\tExonID\tSWscore\tCTCF\tPromoterLike\tProximalEnhancerLike\tDistalEnhancerLike\tDNaseH3K4me3\tintOGEN\tOncoKB\tDNArepair\tmutType\n" #>$NNN\.header.txt

#echo -en "chr\tstart\tend\tCaseID\tref\talt\ttotDepth_B\trefDepth_B\taltDepth_B\ttotDepth_T\trefDepth_T\taltDepth_T\tAFreq_T\trefDepth_R\taltDepth_R\tGeneName\tAAchange\tSnpID\tsnpAF\tCosmicID\tnCosmic\tEffect\tImpact\tFunctionClass\tAAlength\tBioType\tGeneCoding\tTranscriptID\tExonID\tSWscore\tBcell_CTCF\tBcell_H3K27ac\tBcell_SE\tccREs_CTCF\tccREs_PromoterLike\tccREs_ProximalEnhancerLike\tccREs_DistalEnhancerLike\tccREs_DNaseH3K4me3\tintOGEN\tOncoKB\tDNArepair\tmutType\n"

#echo -en "chr\tstart\tend\tCaseID\tref\talt\ttotDepth_B\ttotDepth_T\tAFreq_B\tAFreq_T\tGeneName\tAAchange\tSnpID\tsnpAF\tCosmicID\tnCosmic\tEffect\tImpact\tFunctionClass\tAAlength\tBioType\tGeneCoding\tTranscriptID\tExonID\tClinVar_Significance\tnonCancer_AF\tfold\tfile\tGorS\tmutType\n"

#head -n1 header.txt | awk '{print $0"\tmutType"}'

cat $NNN | awk -F"\t" '{if(length($4)==1 && length($5)==1){print $0"\tSNP"}else if(length($4)>1 && length($5)==1){print $0"\tDel"}else if(length($4)==1 && length($5)>1){print $0"\tIns"}else{print $0"\tINDEL"}}'
