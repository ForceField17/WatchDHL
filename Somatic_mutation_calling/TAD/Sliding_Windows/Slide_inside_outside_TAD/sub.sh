cut -f 1,2,3,4 ../../get_clean_TAD/Final_direction_TAD_boundaries.bed | grep -v 105535000 > Input_peaks.bed # IGH region was removed because of lacking TAD information

sh Density_nege.sh DHL DDD | sed s/,/\\t/g > mutDens_DHL_nege.txt &
sh Density.sh DHL AAA | sed s/,/\\t/g > mutDens_DHL.txt &
sh Density_nege.sh FL CCC | sed s/,/\\t/g > mutDens_FL_nege.txt &
sh Density.sh FL BBB | sed s/,/\\t/g > mutDens_FL.txt &
wait
