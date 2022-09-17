bedtools intersect -a Slide_Both_H3K4me3_and_H3K27ac/Input_peaks.bed -b ComB.bed > Both_H3K4me3_and_H3K27ac_B
bedtools intersect -a Slide_H3K27ac/Input_peaks.bed -b ComB.bed > H3K27ac_B
bedtools intersect -a Slide_H3K4me1/Input_peaks.bed -b ComB.bed > H3K4me1_B
bedtools intersect -a Slide_H3K4me3/Input_peaks.bed -b ComB.bed > H3K4me3_B
bedtools intersect -a Slide_Only_H3K27ac/Input_peaks.bed -b ComB.bed > Only_H3K27ac_B
bedtools intersect -a Slide_Only_H3K4me3/Input_peaks.bed -b ComB.bed > Only_H3K4me3_B



cp Both_H3K4me3_and_H3K27ac_B Slide_Both_H3K4me3_and_H3K27ac/Input_peaks.bed
cp H3K27ac_B Slide_H3K27ac/Input_peaks.bed
cp H3K4me1_B Slide_H3K4me1/Input_peaks.bed
cp H3K4me3_B Slide_H3K4me3/Input_peaks.bed
cp Only_H3K27ac_B Slide_Only_H3K27ac/Input_peaks.bed
cp Only_H3K4me3_B Slide_Only_H3K4me3/Input_peaks.bed

