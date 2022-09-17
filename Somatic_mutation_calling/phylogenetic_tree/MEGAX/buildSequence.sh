echo -en ">P8.germline\n"
cat ../../step4.*csv | sort | uniq | awk '($4=="P8" && $47>=15 && $48>=15 && $49>=15 && $50>=15)' | cut -f 4,5,6,7,8,47,48,49,50,52,53,54,55,108 | awk '{if($10>=10){print $5"\t1\t"$14}else{print $4"\t0\t"$14}}' | awk '{if($3=="Del" && $2==1){print "-"}else if($3=="Del" && $2==0){print "A"}else if($3=="Ins" && $2==1){print "A"}else if($3=="Ins" && $2==0){print "-"}else{print $1}}' | tr "\n" "," | sed s/\,//g
echo -en "\n>P8.a.FL\n" 
cat ../../step4.*csv | sort | uniq | awk '($4=="P8" && $47>=15 && $48>=15 && $49>=15 && $50>=15)' | cut -f 4,5,6,7,8,47,48,49,50,52,53,54,55,108 | awk '{if($11>=10){print $5"\t1\t"$14}else{print $4"\t0\t"$14}}' | awk '{if($3=="Del" && $2==1){print "-"}else if($3=="Del" && $2==0){print "A"}else if($3=="Ins" && $2==1){print "A"}else if($3=="Ins" && $2==0){print "-"}else{print $1}}' | tr "\n" "," | sed s/\,//g
echo -en "\n>P8.b.FL\n"
cat ../../step4.*csv | sort | uniq | awk '($4=="P8" && $47>=15 && $48>=15 && $49>=15 && $50>=15)' | cut -f 4,5,6,7,8,47,48,49,50,52,53,54,55,108 | awk '{if($12>=20){print $5"\t1\t"$14}else{print $4"\t0\t"$14}}' | awk '{if($3=="Del" && $2==1){print "-"}else if($3=="Del" && $2==0){print "A"}else if($3=="Ins" && $2==1){print "A"}else if($3=="Ins" && $2==0){print "-"}else{print $1}}' | tr "\n" "," | sed s/\,//g
echo -en "\n>P8.DHL\n"
cat ../../step4.*csv | sort | uniq | awk '($4=="P8" && $47>=15 && $48>=15 && $49>=15 && $50>=15)' | cut -f 4,5,6,7,8,47,48,49,50,52,53,54,55,108 | awk '{if($13>=20){print $5"\t1\t"$14}else{print $4"\t0\t"$14}}' | awk '{if($3=="Del" && $2==1){print "-"}else if($3=="Del" && $2==0){print "A"}else if($3=="Ins" && $2==1){print "A"}else if($3=="Ins" && $2==0){print "-"}else{print $1}}' | tr "\n" "," | sed s/\,//g
echo -en "\n"
