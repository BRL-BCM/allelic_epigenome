use strict;
use warnings;
use Bio::DB::Sam;

my $bamFile = shift(@ARGV);
my $refGenome = shift(@ARGV);
my $roisFile = shift(@ARGV);

my $sam = Bio::DB::Sam->new(-bam=>$bamFile,
			    -fasta=>$refGenome,
			    -autoindex=>1);

open(ROIS,$roisFile);
print "Chromosome","\t",
    "Start","\t",
    "End","\t",
    "CpG_Positions","\t",
    "Het_Positions","\t",
    "Patterns","\n";

while(my $roi = <ROIS>){
    chomp($roi);
    if($roi =~ /^\s*$/){
	next;
    }
    my ($roiChr,$roiStart,$roiEnd,$cpgPos,$cpgStrands,$hetPos,$als1,$als2) = split("\t",$roi);
    
    if($cpgPos eq "."){next;}
    my @cpgPos = split(",",$cpgPos);
    my @cpgStrands = split(",",$cpgStrands);
    if($cpgStrands[0] eq "-"){
	$cpgStrands = "+".$cpgStrands;
	my $temp = $cpgPos[0]-1;
	$cpgPos = $temp.$cpgPos;
	unshift(@cpgPos,$temp);
	unshift(@cpgStrands,"+");
    }
    if($cpgStrands[-1] eq "+"){
	$cpgStrands = $cpgStrands."-";
	my $temp = $cpgPos[-1]+1;
	$cpgPos = $cpgPos.$temp;
	push(@cpgPos,$temp);
	push(@cpgStrands,"-");
    }
    my $nCpgs = scalar(@cpgPos);
    my @cpgPosOut = @cpgPos[grep($_%2==0, 0..$#cpgPos)]; #get only the positions in the forward strand
    my $cpgPosOut = join(",",@cpgPosOut);
    my $nCpgsOut = scalar(@cpgPosOut);
    my $curCpgOut = 0;
    
    my @hetPos = split(",",$hetPos);
    my $nHets = scalar(@hetPos);
    my @als1 = split(",",$als1);
    my @als2 = split(",",$als2);
    if($hetPos eq "."){$nHets = 0;}
    
    my %readsHets;
    my %readsCpgs;
    my $curHet = 0;
    my $curCpg = 0;
    my $prevPos = 0;

    my $getRoiReads = sub {
	my ($chr,$pos,$pileups) = @_[0..2];
	
	#if($prevPos != 0 && $prevPos != $pos-1){
	#    print STDERR "Warning: Skipped a position!\n";
	#}	
	#$prevPos = $pos;

	if($curHet < @hetPos && $nHets>0){
	    if($hetPos[$curHet] < $pos){
		$curHet++;
	    }
	}
	if($curCpg < @cpgPos){
	    if($cpgPos[$curCpg] < $pos){
		if($cpgStrands[$curCpg] eq "-"){
		    $curCpgOut++;
		}
		$curCpg++;
	    }
	}

	if($nHets>0 && $curHet < @hetPos){
	    if($pos == $hetPos[$curHet]){
		my $al1 = $als1[$curHet];
		my $al2 = $als2[$curHet];
		my $noNegAs = 0;
		my $noPosTs = 0;
		my $orAl1 = $al1;
		my $orAl2 = $al2;
		if(($al1 =~ /A/i and $al2 =~ /C/i) or ($al2 =~ /A/i and $al1 =~ /C/i)){
		    $al1 = "A";
		    $al2 = "[CT]";
		    $orAl1 = "A";
		    $orAl2 = "C";
		}
		elsif(($al1 =~ /T/i and $al2 =~ /G/i) or ($al2 =~ /T/i and $al1 =~ /G/i)){
		    $al1 = "T";
		    $al2 = "[GA]";
		    $orAl1 = "T";
		    $orAl2 = "G";
		}
		elsif(($al1 =~ /C/i and $al2 =~ /G/i) or ($al2 =~ /C/i and $al1 =~ /G/i)){
		    $al1 = "[CT]";
		    $al2 = "[GA]";
		    $orAl1 = "C";
		    $orAl2 = "G";
		}
		elsif(($al1 =~ /A/i and $al2 =~ /G/i) or ($al2 =~ /A/i and $al1 =~ /G/i)){
		    $noNegAs = 1;
		}
		elsif(($al1 =~ /C/i and $al2 =~ /T/i) or ($al2 =~ /C/i and $al1 =~ /T/i)){
		    $noPosTs = 1;
		}
		elsif(($al1 =~ /C/i and $al2 =~ /C/i)){
		    $al1 = "[CT]";
		    $al2 = "[CT]";
		    $orAl1 = "C";
		    $orAl2 = "C";
		}
		elsif(($al1 =~ /G/i and $al2 =~ /G/i)){
		    $al1 = "[GA]";
		    $al2 = "[GA]";
		    $orAl1 = "G";
		    $orAl2 = "G";
		}
		
		for my $pileup (@$pileups) {
		    if($pileup->indel or $pileup->is_refskip){
			next; # Ignore reads with indels in this position
		    }
		    my $alignment = $pileup->alignment;
		    my $qscore = $alignment->qscore->[$pileup->qpos];
		    if($qscore < 20){
			next; # Ignore reads with base quality less than 20 in this position
		    }
		    my $qbase  = substr($alignment->qseq,$pileup->qpos,1);
		    my $strand = $alignment->strand;
		    
		    if($noPosTs and (($qbase =~ /T/i) && $strand == 1)){
			next; # Ignore reads with T basecall from forward when genotype is C/T
		    }
		    elsif($noNegAs and  (($qbase =~ /A/i) && $strand == -1)){
			next; # Ignore reads with basecall A from reverse strand when genotype is A/G
		    }
		    my $readName = $alignment->qname;
		    unless(exists($readsHets{$readName})){
			$readsHets{$readName} = [("*") x $nHets]; 
		    } 
		    if($qbase =~ /$al1/i){
			$readsHets{$readName}->[$curHet] = $orAl1;
		    }
		    elsif($qbase =~ /$al2/i){
			$readsHets{$readName}->[$curHet] = $orAl2;
		    }
		}
		$curHet++;
	    }
	}
	if($curCpg < @cpgPos){
	    if($pos == $cpgPos[$curCpg]){
		my $cStrand = $cpgStrands[$curCpg];
		for my $pileup (@$pileups){
		    if($pileup->indel or $pileup->is_refskip){
			next;
		    }
		    my $alignment = $pileup->alignment;
		    my $qscore = $alignment->qscore->[$pileup->qpos];
		    if($qscore < 20){
			next;
		    }
		    my $qbase  = substr($alignment->qseq,$pileup->qpos,1);
		    my $strand = $alignment->strand;
		    my $readName = $alignment->qname;
		    unless(exists($readsCpgs{$readName})){
			$readsCpgs{$readName} = [(0) x $nCpgsOut];
		    } 
		    if(($qbase=~ /C/i && $strand == 1 && $cStrand eq "+")){
			$readsCpgs{$readName}->[$curCpgOut] = 2;
		    }
		    elsif(($qbase=~ /G/i && $strand == -1 && $cStrand eq "-")){
			$readsCpgs{$readName}->[$curCpgOut] = 2;
		    }
		    elsif(($qbase=~ /T/i && $strand == 1 && $cStrand eq "+")){
			$readsCpgs{$readName}->[$curCpgOut] = 1;
		    }
		    elsif(($qbase=~ /A/i && $strand == -1 && $cStrand eq "-")){
			$readsCpgs{$readName}->[$curCpgOut] = 1;
		    }
		}
		$curCpg++;
		if($cStrand eq "-"){$curCpgOut++;}
	    }
	}
    };
    $sam->fast_pileup("$roiChr:$roiStart-$roiEnd",$getRoiReads);
    my %patterns;
    foreach my $readName (keys(%readsCpgs)){
	my $pattern = join("",@{$readsCpgs{$readName}});
	if(exists($readsHets{$readName})){
	    $pattern = $pattern.":".join("",@{$readsHets{$readName}});
	}
	else{
	    $pattern = $pattern.":".join("",("*")x$nHets);
	}
	if(!exists($patterns{$pattern})){
	    $patterns{$pattern} = 1;
	}
	else{
	    $patterns{$pattern}++;
	}
    }
    
    print $roiChr,"\t",$roiStart,"\t",$roiEnd,"\t",$cpgPosOut,"\t",$hetPos,"\t";
    my $isFirst = 1;
    foreach my $pattern (sort(keys(%patterns))){
	if($isFirst){
	    print $pattern,":",$patterns{$pattern};
	    $isFirst = 0;
	}
	else{
	    print ",",$pattern,":",$patterns{$pattern};
	}
    }
    if($isFirst == 1){print ".";}
    print "\n";
}
