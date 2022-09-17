#!/usr/bin/perl

#please run this job in a new clear fold
die "usage: perl checking.pl <savi_table> \n" if @ARGV!= 1;

#`sh grep1.sh | sed s/SE/\\tSE/g | sed s/' '//g >top20SE.txt`;

$target = $ARGV[0];
`sh grep1.sh >$target`;
open FILE, $target;
print "gene\tproteinLength\tgeneCoding\tchr\tOncoKB\tpublicationReported\tNumRecordsAll\tNumR_FL\tNumR_DHL\tNumPatientsAll\tNumP_FL\tNumP_DHL\tCosmicRecord\tImpact\tClinVar\tMutatedFL\tMutatedDHL\n";
$i=0;
while($line = <FILE>){
	chomp($line);
	#$content[$i] = $line;
	@temp=split('\s+',$line);
	$nMut = $temp[1];
	$H3 = $temp[0];
	$SE1 = $temp[2];
        $SE2 = $temp[3];
        $SE3 = $temp[4];
	$SE4 = $temp[5];
	$SE5 = $temp[6];
	$xx1 = `sh grep2NR.sh "$H3" temp`;
	chomp($xx1);
        $xx2 = `sh grep2NR.sh "$H3" step4.FL.csv`;
        chomp($xx2);
        $xx3 = `sh grep2NR.sh "$H3" step4.DHL.csv`;
        chomp($xx3);

	$xx4 = `sh grep2NP.sh "$H3" temp`;
        chomp($xx4);
        $xx5 = `sh grep2NP.sh "$H3" step4.FL.csv`;
        chomp($xx5);
        $xx6 = `sh grep2NP.sh "$H3" step4.DHL.csv`;
        chomp($xx6);

	$genes = `sh grep3.sh "$H3" temp`;
        chomp($genes);
        $genes =~ s/\,$//g;

        $hancer = `sh grep4.sh "$H3" temp`;
        chomp($hancer);
        $hancer =~ s/\,$//g;

	$histon = `sh grep5.sh "$H3" temp`;
        chomp($histon);
        $histon =~ s/\,$//g;

	$xx7 = `sh grep6.sh "$H3" step4.FL.csv`;
        chomp($xx7);
        $xx8 = `sh grep6.sh "$H3" step4.DHL.csv`;
        chomp($xx8);

	print "$H3\t$SE1\t$SE2\t$SE3\t$SE4\t$SE5\t$xx1\t$xx2\t$xx3\t$xx4\t$xx5\t$xx6\t$genes\t$hancer\t$histon\t$xx7\t$xx8\n";	
}
close FILE;


