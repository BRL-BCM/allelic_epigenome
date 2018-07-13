use strict;
use warnings;

my $vcfFile = shift(@ARGV);
my $aiTable = shift(@ARGV);

open(VCF, $vcfFile);

open(AIT, $aiTable);

my $header = <AIT>;
chomp($header);
$header = $header."\trsID\tPhasing.Block.Start\tPhasing.Orientation";
print $header,"\n";

while(my $aitLine = <AIT>){
    my $vcfLine = <VCF>;
    while($vcfLine =~ /^#/ || $vcfLine =~ /1\/1/){
	$vcfLine = <VCF>;
    }
    chomp($vcfLine);
    chomp($aitLine);
    my @vcfLine = split("\t",$vcfLine);
    my @aitLine = split("\t",$aitLine);
    
    if(!($vcfLine[0] eq $aitLine[0] && $vcfLine[1] eq $aitLine[2])){
	print STDERR "Vcf file line does not match the allelic imbalance table line\n";
	print STDERR "VCF line: $vcfLine\n";
	print STDERR "AIT line: $aitLine\n";
	#$aitLine = <AIT>;
	#@aitLine = split("\t",$aitLine);
	last;
    }

    my $alleleFlip = 1;
    if(!($aitLine[9] eq $vcfLine[3] && $aitLine[10] eq $vcfLine[4])){
	if($aitLine[10] eq $vcfLine[3] && $aitLine[9] eq $vcfLine[4]){
	    $alleleFlip = -1;
	}
	else{
	    print STDERR "Alleles in VCF do not match alleles in allelic imbalance table\n";
	    last;
	}
    }
	

    my @vcfInfo = split(":",$vcfLine[9]);
    my ($phasing1,$phasing2) = split(",",$vcfInfo[2]);
    my ($phasingBlock,$orientation) = split("-",$phasing1);

    if($orientation == 2){
	$orientation = "-1";
    }
    $orientation *= $alleleFlip;
    print $aitLine,"\t$vcfLine[2]\t$aitLine[0]:$phasingBlock\t$orientation\n";
}

close(VCF);
close(AIT);
