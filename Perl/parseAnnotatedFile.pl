use strict;
use warnings;

print "#chromosome","\t",
    "start","\t",
    "end","\t",
    "allele1.methylated","\t",
    "allele1.unmethylated","\t",
    "allele2.methylated","\t",
    "allele2.unmethylated","\t",
    "total.number.of.reads","\t",
    "maximum.number.of.cpgs","\t",
    "number.of.heterozygous.samples","\t",
    "allele1","\t",
    "allele2","\t",
    "rsID","\t",
    "allele1.frequency","\t",
    "allele2.frequency","\t",
    "ancestral.allele","\t",
    "1kg.pilot.mask","\t",
    "1kg.strict.mask","\t",
    "is.reference.allele.present","\n";

my $lineNumber = 0;
while(my $line = <>){
    $lineNumber++;
    chomp($line);
    my @fields = split("\t",$line);
    my $al1 = $fields[10];
    my $al2 = $fields[11];
    my $refAl = $fields[13];
    my $altAls = $fields[14];
    my $alFreqs = $fields[15];

    my $rsId = $fields[12];
    if($rsId =~ /,/){
	#print STDERR "Line $lineNumber has multiple rsIDs\n";
	#print STDERR $line,"\n";
	my @temp = split(",",$rsId);
	$rsId = $temp[0];
    }

    my @possibleAls;
    push(@possibleAls,$refAl);
    if($altAls =~ /,/){
	my @altAls = split(",",$altAls);
	push(@possibleAls,@altAls);
    }
    else{
	push(@possibleAls,$altAls);
    }

    my @alFreqs;
    if($alFreqs =~ /,/){
	@alFreqs = split(",",$alFreqs);
    }
    else{
	@alFreqs = (".") x scalar(@possibleAls);
    }

    if(scalar(@possibleAls) != scalar(@alFreqs)){
	print STDERR "Line $lineNumber has strange set of alleles/frequencies\n";
    }

    my $al1Freq = ".";
    my $al2Freq = ".";
    for(my $i = 0; $i < @possibleAls; $i++){
	if($al1 eq $possibleAls[$i]){
	    $al1Freq = $alFreqs[$i];
	}
	elsif($al2 eq $possibleAls[$i]){
	    $al2Freq = $alFreqs[$i];
	}
    }
    
    if($al1 eq $refAl){
	print join("\t",@fields[0..11]),"\t",$rsId,"\t",$al1Freq,"\t",$al2Freq,"\t",join("\t",@fields[16..18]),"\t",1,"\n";
    }
    elsif($al2 eq $refAl){
	print join("\t",@fields[0..2]),"\t",$fields[5],"\t",$fields[6],"\t",$fields[3],"\t",$fields[4],"\t",join("\t",@fields[7..9]),"\t",$al2,"\t",$al1,"\t",$rsId,"\t",$al2Freq,"\t",$al1Freq,"\t",join("\t",@fields[16..18]),"\t",1,"\n";
    }
    else{
	print join("\t",@fields[0..11]),"\t",$rsId,"\t",$al1Freq,"\t",$al2Freq,"\t",join("\t",@fields[16..18]),"\t",0,"\n";
    }
}
    
    
