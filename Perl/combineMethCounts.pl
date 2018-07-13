use strict;
use warnings;
use POSIX;

my $fileList = shift(@ARGV);
my $minCountsPerAllele = shift(@ARGV);
open(LIST,$fileList);

my @inputFileHandles;
while(my $inFileName = <LIST>){
    chomp($inFileName);
    local *INFILE;
    open(INFILE, $inFileName) || die;
    <INFILE>;
    push(@inputFileHandles, *INFILE);
}

print "#Chromosome","\t",
    "Start","\t",
    "End","\t",
    "Allele1/methylated","\t",
    "Allele1/unmethylated","\t",
    "Allele2/methylated","\t",
    "Allele2/unmethylated","\t",
    "Total.Number.Of.Reads","\t",
    "Maximum.Number.Of.CpGs","\t",
    "Number.Of.Heterozygous.Samples","\t",
    "Allele1","\t",
    "Allele2","\n";

my $fh1 = $inputFileHandles[0];
while(my $lineFile1 = <$fh1>){
    chomp($lineFile1);
    my @lines;
    my @fields = split("\t",$lineFile1);
    my $chr = $fields[0];
    my $start = $fields[1];
    my $end = $fields[2];
    my %totalCounts;
    my $totalReads = 0;
    my @nCpGs;
    my $nHetSamples = 0;
    
    my $al1 = $fields[9];
    my $al2 = $fields[10];
    if(!($al1 eq $al2)){
	$nHetSamples++;
    }
    if(!exists($totalCounts{$al1})){
	$totalCounts{$al1} = [0,0];
    }
    if(!exists($totalCounts{$al2})){
	$totalCounts{$al2} = [0,0];
    }
    
    $totalCounts{$al1}->[0] += $fields[3];
    $totalCounts{$al1}->[1] += $fields[4];
    $totalCounts{$al2}->[0] += $fields[5];
    $totalCounts{$al2}->[1] += $fields[6];
    
    $totalReads += $fields[7];
    push(@nCpGs,$fields[8]);
    
    for(my $i = 1; $i < @inputFileHandles; $i++){
	my $fh = $inputFileHandles[$i];
	my $line = <$fh>;
	chomp($line);
	@fields = split("\t",$line);
	$al1 = $fields[9];
	$al2 = $fields[10];
	if(!($al1 eq $al2)){
	    $nHetSamples++;
	}
	if(!exists($totalCounts{$al1})){
	    $totalCounts{$al1} = [0,0];
	}
	if(!exists($totalCounts{$al2})){
	    $totalCounts{$al2} = [0,0];
	}
	
	$totalCounts{$al1}->[0] += $fields[3];
	$totalCounts{$al1}->[1] += $fields[4];
	$totalCounts{$al2}->[0] += $fields[5];
	$totalCounts{$al2}->[1] += $fields[6];
	
	$totalReads += $fields[7];
	push(@nCpGs,$fields[8]);
    }

    my @alleles = keys(%totalCounts);

    if(scalar(@alleles) != 2){
	next;
    }
    
    $al1 = $alleles[0];
    $al2 = $alleles[1];
    @nCpGs = sort(@nCpGs);
    my $nSamples = scalar(@nCpGs);
    my $medianNCpGs = $nCpGs[$nSamples-1];

    if($totalCounts{$al1}->[0]+$totalCounts{$al1}->[1] >= $minCountsPerAllele && $totalCounts{$al2}->[0]+$totalCounts{$al2}->[1] >= $minCountsPerAllele){
	print 
	    $chr,"\t",
	    $start,"\t",
	    $end,"\t",
	    $totalCounts{$al1}->[0],"\t",
	    $totalCounts{$al1}->[1],"\t",
	    $totalCounts{$al2}->[0],"\t",
	    $totalCounts{$al2}->[1],"\t",
	    $totalReads,"\t",
	    $medianNCpGs,"\t",
	    $nHetSamples,"\t",
	    $al1,"\t",
	    $al2,"\n";
    }
}

