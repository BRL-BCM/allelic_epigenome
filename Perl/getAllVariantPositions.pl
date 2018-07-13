use strict;
use warnings;

my $vcfFileList = shift(@ARGV);
my $nVcfFiles = shift(@ARGV);

my %allVars;
open(VCFLIST,$vcfFileList);
my $vcfFileNumber = 0;
my @vcfTags;
while(my $line = <VCFLIST>){
    chomp($line);
    my ($vcfFile,$vcfTag) = split(" ", $line);
    push(@vcfTags,$vcfTag);
    `zcat $vcfFile > currentVcf.vcf`;
    open(VCF,"currentVcf.vcf");
    while(my $line = <VCF>){
	if($line =~ /^#/){
	    next;
	}
	chomp($line);
	my @fields = split("\t",$line);
	my $chr = $fields[0];
	my $pos = $fields[1];
	my $refAl = $fields[3];
	my @altAl;
	if($fields[4] =~ /,/){
	    @altAl = split(",",$fields[4]);
	}
	else{
	    push(@altAl,$fields[4]);
	}
	
	my $otherInfo = $fields[9];
	my @otherInfo = split(":",$otherInfo);
	my $gt = $otherInfo[0];
	if($gt eq "0/0" || !($gt =~ /^\d\/\d$/)){
	    next;
	}
	$gt =~ s/0/$refAl/g;
	for(my $i = 1; $i <= @altAl; $i++){
	    $gt =~ s/$i/$altAl[$i-1]/g;
	}
	my ($al1,$al2) = split("/",$gt);
	
	my $key = $chr."\t".($pos-1)."\t".$pos;
	if(!exists($allVars{$key})){
	    my @temp = ("$refAl\t$refAl") x $nVcfFiles;
	    $allVars{$key} = \@temp;
	}
	
	$allVars{$key}->[$vcfFileNumber] = "$al1\t$al2";
    }
    print STDERR "Done with $vcfTag.\n";
    `rm currentVcf.vcf`;
    $vcfFileNumber++;
}
close(VCFLIST);

for(my $i = 0; $i < @vcfTags; $i++){
    open(ALLVARS,">$vcfTags[$i].allVariantPositions.txt");
    foreach my $var (keys(%allVars)){
	print ALLVARS $var,"\t",$allVars{$var}->[$i],"\n";
    }
    close(ALLVARS);
}
