#P1 BCL2
perl trans.pl TACTACTAGCCGGAAGGAGG > pattern_P1
#P1 MYC
perl trans.pl ACCCACTTTAGCCCGCCTGT >> pattern_P1

zcat ../P1_NL.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P1 > reads1.P1_NL &
zcat ../P1_NL.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P1 > reads2.P1_NL &
zcat ../P1_LG.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P1 > reads1.P1_FL &
zcat ../P1_LG.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P1 > reads2.P1_FL &
zcat ../P1_HG.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P1 > reads1.P1_DHL &
zcat ../P1_HG.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P1 > reads2.P1_DHL &
zcat ../P1.ffpe_HG.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P1 > reads1.P1_DHL_ffpe &
zcat ../P1.ffpe_HG.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P1 > reads2.P1_DHL_ffpe &


#P2 BCL2
perl trans.pl GCGACTAGTTAGCACCACTG > pattern_P2
#P2 MYC
perl trans.pl TCATCTTCTCCTATTCTGCC >> pattern_P2

zcat ../P2_NL.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P2 > reads1.P2_NL &
zcat ../P2_NL.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P2 > reads2.P2_NL &
zcat ../P2_LG.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P2 > reads1.P2_FL &
zcat ../P2_LG.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P2 > reads2.P2_FL &
zcat ../P2_HG.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P2 > reads1.P2_DHL &
zcat ../P2_HG.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P2 > reads2.P2_DHL &

wait


#P3 BCL2
perl trans.pl GGTCCTGGACCGAGCTGTTC > pattern_P3
#P3 MYC
perl trans.pl TAGGGCGATGCCGTCTCCGG >> pattern_P3

zcat ../P3_LG.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P3 > reads1.P3_FL &
zcat ../P3_LG.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P3 > reads2.P3_FL &
zcat ../P3_HG.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P3 > reads1.P3_DHL &
zcat ../P3_HG.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P3 > reads2.P3_DHL &


#P4 BCL2
perl trans.pl ACCAGCTCTTTTTCAGGAA > pattern_P4
#P4 MYC
perl trans.pl TTTGCAAAAGATTTGGTTC >> pattern_P4

zcat ../P4_LG.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P4 > reads1.P4_FL &
zcat ../P4_LG.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P4 > reads2.P4_FL &
zcat ../P4_HG.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P4 > reads1.P4_DHL &
zcat ../P4_HG.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P4 > reads2.P4_DHL &


#P5 BCL2
perl trans.pl GCCGCCCCGAAGGAGGGC > pattern_P5
#P5 MYC
perl trans.pl ATAGCAAGATAGGTTTGT >> pattern_P5

zcat ../P5_LG.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P5 > reads1.P5_FL &
zcat ../P5_LG.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P5 > reads2.P5_FL &
zcat ../P5_HG.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P5 > reads1.P5_DHL &
zcat ../P5_HG.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P5 > reads2.P5_DHL &


#P6 BCL2
perl trans.pl ACGTGGGCTTTGCAGCAT > pattern_P6
#P6 MYC
perl trans.pl CCAAACCAACTGGGGCTG >> pattern_P6

zcat ../P6_NL.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P6 > reads1.P6_NL &
zcat ../P6_NL.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P6 > reads2.P6_NL &
zcat ../P6_LG.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P6 > reads1.P6_FL &
zcat ../P6_LG.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P6 > reads2.P6_FL &
zcat ../P6_HG.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P6 > reads1.P6_DHL &
zcat ../P6_HG.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P6 > reads2.P6_DHL &


#P7 BCL2
perl trans.pl AGAGGTGGGCTTCATACC > pattern_P7
#P7 MYC
perl trans.pl AAAAGCGGGAGCAAACAA >> pattern_P7

zcat ../P7_NL.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P7 > reads1.P7_NL &
zcat ../P7_NL.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P7 > reads2.P7_NL &
zcat ../P7_LG.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P7 > reads1.P7_FL &
zcat ../P7_LG.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P7 > reads2.P7_FL &
zcat ../P7_HG.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P7 > reads1.P7_DHL &
zcat ../P7_HG.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P7 > reads2.P7_DHL &


#P8 BCL2
perl trans.pl CACCATCTTAAGCACC > pattern_P8
#P8 MYC
perl trans.pl GTTCGCTGGTGACCCCAG >> pattern_P8

zcat ../P8_NL.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P8 > reads1.P8_NL &
zcat ../P8_NL.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P8 > reads2.P8_NL &
zcat ../P8.a_LG.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P8 > reads1.P8_FL.a &
zcat ../P8.a_LG.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P8 > reads2.P8_FL.a &
zcat ../P8.b_LG.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P8 > reads1.P8_FL.b &
zcat ../P8.b_LG.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P8 > reads2.P8_FL.b &
zcat ../P8_HG.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P8 > reads1.P8_DHL &
zcat ../P8_HG.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P8 > reads2.P8_DHL &


#P9 BCL2
perl trans.pl AAACATTAAGCACCACTG > pattern_P9
#P9 MYC
perl trans.pl GTTGGGAGCGTAGTTT >> pattern_P9

zcat ../P9_NL.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P9 > reads1.P9_NL &
zcat ../P9_NL.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P9 > reads2.P9_NL &
zcat ../P9_LG.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P9 > reads1.P9_FL &
zcat ../P9_LG.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P9 > reads2.P9_FL &
zcat ../P9_HG.R1.fastq.gz | LC_ALL=C fgrep -f pattern_P9 > reads1.P9_DHL &
zcat ../P9_HG.R2.fastq.gz | LC_ALL=C fgrep -f pattern_P9 > reads2.P9_DHL &

wait

#integrate and unique all splite reads
mkdir Unique_Reads
cat reads1.P1_NL reads2.P1_NL | sort | uniq > Unique_Reads/ReadsPool.P1_NL
cat reads1.P2_NL reads2.P2_NL | sort | uniq > Unique_Reads/ReadsPool.P2_NL
cat reads1.P6_NL reads2.P6_NL | sort | uniq > Unique_Reads/ReadsPool.P6_NL
cat reads1.P7_NL reads2.P7_NL | sort | uniq > Unique_Reads/ReadsPool.P7_NL
cat reads1.P8_NL reads2.P8_NL | sort | uniq > Unique_Reads/ReadsPool.P8_NL
cat reads1.P9_NL reads2.P9_NL | sort | uniq > Unique_Reads/ReadsPool.P9_NL

cat reads1.P1_FL   reads2.P1_FL   | sort | uniq > Unique_Reads/ReadsPool.P1_FL
cat reads1.P2_FL   reads2.P2_FL   | sort | uniq > Unique_Reads/ReadsPool.P2_FL
cat reads1.P3_FL   reads2.P3_FL   | sort | uniq > Unique_Reads/ReadsPool.P3_FL
cat reads1.P4_FL   reads2.P4_FL   | sort | uniq > Unique_Reads/ReadsPool.P4_FL
cat reads1.P5_FL   reads2.P5_FL   | sort | uniq > Unique_Reads/ReadsPool.P5_FL
cat reads1.P6_FL   reads2.P6_FL   | sort | uniq > Unique_Reads/ReadsPool.P6_FL
cat reads1.P7_FL   reads2.P7_FL   | sort | uniq > Unique_Reads/ReadsPool.P7_FL
cat reads1.P8_FL.a reads2.P8_FL.a | sort | uniq > Unique_Reads/ReadsPool.P8_FL_a
cat reads1.P8_FL.b reads2.P8_FL.b | sort | uniq > Unique_Reads/ReadsPool.P8_FL_b
cat reads1.P9_FL   reads2.P9_FL   | sort | uniq > Unique_Reads/ReadsPool.P9_FL

cat reads1.P1_DHL   reads2.P1_DHL   | sort | uniq > Unique_Reads/ReadsPool.P1_DHL
cat reads1.P2_DHL   reads2.P2_DHL   | sort | uniq > Unique_Reads/ReadsPool.P2_DHL
cat reads1.P3_DHL   reads2.P3_DHL   | sort | uniq > Unique_Reads/ReadsPool.P3_DHL
cat reads1.P4_DHL   reads2.P4_DHL   | sort | uniq > Unique_Reads/ReadsPool.P4_DHL
cat reads1.P5_DHL   reads2.P5_DHL   | sort | uniq > Unique_Reads/ReadsPool.P5_DHL
cat reads1.P6_DHL   reads2.P6_DHL   | sort | uniq > Unique_Reads/ReadsPool.P6_DHL
cat reads1.P7_DHL   reads2.P7_DHL   | sort | uniq > Unique_Reads/ReadsPool.P7_DHL
cat reads1.P8_DHL   reads2.P8_DHL   | sort | uniq > Unique_Reads/ReadsPool.P8_DHL
cat reads1.P9_DHL   reads2.P9_DHL   | sort | uniq > Unique_Reads/ReadsPool.P9_DHL




