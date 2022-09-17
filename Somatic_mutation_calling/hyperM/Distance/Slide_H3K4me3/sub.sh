sh Density_nege.sh DHL DDD | sed s/,/\\t/g > mutDens_DHL_nege.txt &
sh Density.sh DHL AAA | sed s/,/\\t/g > mutDens_DHL.txt &
sh Density_nege.sh FL CCC | sed s/,/\\t/g > mutDens_FL_nege.txt &
sh Density.sh FL BBB | sed s/,/\\t/g > mutDens_FL.txt &
wait
