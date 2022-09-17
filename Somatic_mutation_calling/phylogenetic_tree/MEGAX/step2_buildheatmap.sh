#mkdir heatmap
head -n1  ../../step4.FL.csv | cut -f 4,5,6,7,8,47,48,49,50,52,53,54,55,108 | awk '{print $0"\tgermline\tFL.a\tFL.b\tDHL"}' 
cat ../../step4.*csv | sort | uniq | awk '($4=="P8" && $47>=15 && $48>=15 && $49>=15 && $50>=15)' | cut -f 4,5,6,7,8,47,48,49,50,52,53,54,55,108 | awk '{if($10>=10){print $0"\tMut"}else{print $0"\tWT"}}'   | awk '{if($11>=10){print $0"\tMut"}else{print $0"\tWT"}}'  | awk '{if($12>=20){print $0"\tMut"}else{print $0"\tWT"}}'  | awk '{if($13>=20){print $0"\tMut"}else{print $0"\tWT"}}' 
