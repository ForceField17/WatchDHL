bigWigToWig 4DNFILYQ1PAY_compartmentAB.bw compartment.wig
cat compartment.wig | awk '($4>0.01 && $4!="nan")' > AAA
bedtools merge -i AAA | awk '{print $0"\tA"}' > CompartmentA.bed

cat compartment.wig | awk '($4<-0.01 && $4!="nan")' > BBB
bedtools merge -i BBB | awk '{print $0"\tB"}' > CompartmentB.bed

cat CompartmentA.bed CompartmentB.bed > CCC

sort-bed CCC > Final_compartment.bed
rm AAA BBB CCC

