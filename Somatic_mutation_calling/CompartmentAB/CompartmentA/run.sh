bedtools intersect -a Slide_Both_H3K4me3_and_H3K27ac/Input_peaks.bed -b ComA.bed > Both_H3K4me3_and_H3K27ac_A
bedtools intersect -a Slide_H3K27ac/Input_peaks.bed -b ComA.bed > H3K27ac_A
bedtools intersect -a Slide_H3K4me1/Input_peaks.bed -b ComA.bed > H3K4me1_A
bedtools intersect -a Slide_H3K4me3/Input_peaks.bed -b ComA.bed > H3K4me3_A
bedtools intersect -a Slide_Only_H3K27ac/Input_peaks.bed -b ComA.bed > Only_H3K27ac_A
bedtools intersect -a Slide_Only_H3K4me3/Input_peaks.bed -b ComA.bed > Only_H3K4me3_A



cp Both_H3K4me3_and_H3K27ac_A Slide_Both_H3K4me3_and_H3K27ac/Input_peaks.bed
cp H3K27ac_A Slide_H3K27ac/Input_peaks.bed
cp H3K4me1_A Slide_H3K4me1/Input_peaks.bed
cp H3K4me3_A Slide_H3K4me3/Input_peaks.bed
cp Only_H3K27ac_A Slide_Only_H3K27ac/Input_peaks.bed
cp Only_H3K4me3_A Slide_Only_H3K4me3/Input_peaks.bed

