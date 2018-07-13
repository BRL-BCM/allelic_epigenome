# Code related to Allelic Epigenome 

Make sure reference genome files (.fa files) and BAM files are indexed and sorted prior to running the following codes.

**1**. To call allele-specific methylation from a single sample

**1.1**	To generate a file containing positions of homozygous CpGs and heterozygous SNP positions use the getCpgPositions.pl script. The script takes as input the relevant reference genome fasta file (.fa), a variant call file (.vcf) from the sample of interest, and a job name. Returned are two .txt files, one containing the heterozygous SNPs (.hetSnps.txt.gz) and another for homoCpGs (.homoCpGs.txt.gz).  Note: In the hetSnps file, a value of 1 marks heterozygous SNP positions that overlap a CpG position.  

**1.2** To generate an allelic methylation count file use the getAllelicMethCounts.pl script. As input supply the BAM file from your bisulfite-sequencing experiment, the relevant fasta file for your reference genome, and the two previously generated .txt files containing the homozygous CpG positions (from Step 1.1) and heterozygous SNP positions (from Step 1.1).  Sorting input files by chromosome and start position prior to running the script will result in faster run time. Returned is the allelic methylation count .txt file.  

**1.3**	To call allele-specific methylation from the allelic methylation count file use the processMethCounts.R script. Input is the allelic methylation count .txt file (from Step 1.2), the desired name of the output file, and an FDR value cut-off (default 0.1). 

**1.4**	To add phasing information to the allele-specific methylation file (from Step 1.3) use the addPhasingInfo.pl script. Inputs are the variant call file (.vcf) from the sample of interest and allele-specific methylation file (from Step 1.3). 

**2**. Calling SD-ASM in dataset combining reads across multiple samples at positions with 2 total alleles across all donors

**2.1** To get all variant positions across all donors from which reads are combined use the getAllVariantPositions.pl script. Input a .txt file containing two columns separated by a space. The first column should contain the path to the VCF files for each donor and in the second column an identifier for the donor should be created. Returned will be an all variant position file for each donor. 

**2.2** To make an allelic methylation count file across relevant positions for each donor use the getAllelicMethCountsIncludingHomozygousPositions.pl script. Input the bisulfite-sequencing BAM file, the relevant reference genome fasta file, the path to the homozygous CpG file for that donor (from Step 1.1), the all variant position file (from Step 2.1). 

**2.3** To combine allelic methylation count files generated (in Step 2.2) use the combineMethCounts.pl script. The inputs are a .txt file containing the paths to all the allelic methylation count files (from Step 2.2) in separate rows and an integer to specify the minimum number of counts per allele you want to use (e.g. 6 CpG counts used for the paper).

**2.4** To add additional data (such 1000 genomes frequencies, rsIDs, ancestral allele information, etc.) combine the combined meth count file (from Step 2.3) with the annotation.txt by chromosome and position found on http://genboree.org/theCommons/documents/788. The parseAnnotatedFile.pl script will then parse the combined meth count file (with the added annotation data) into a cleaner format as output. Input is the combined meth count file with the annotation data added. 

**2.5** To call allele-specific methylation on the parsed combined methylation count file, use the processCombinedMethCountsFunctions.R script. Relevant functions and their description are commented within the script. Input is the parsed the combined meth count file (from Step 2.4). 

**3.** Calling epialleles from BAM files generated from bisulfite sequencing experiments.

**3.1** To process the BAM file and extract patterns of methylation for windows of 4 consecutive CpGs, utilize the getMethStates.pl script.  Inputs are the BAM file, a fasta file of the relevant reference genome, and a BED file containing the regions of interest. The BED file should contain the chromosome positions in the first column, the coordinate of the start of the window of interest in the second column, the coordinate of the end of the window of interest in the third column, single coordinates of any homozygous CpGs separated by commas in the fourth column (e.g. 12331,10541), the strand of the CpGs separated by commas in the fifth column (e.g +,+), the coordinates of heterozygous SNP positions in the file in the sixth column separated by commas (e.g. 12345,12345), the nucleotides of the first alleles separated by commas should be in the seventh column (e.g. A,C) , and the nucleotides of the second alleles separated by commas should be in the eighth column (e.g. T,G). Returned will be a methylation state file that can be further evaluated using the following R functions.  Note: The VCF files and homozygous CpG files (from Step 1.1) for the sample can be utilized with the mapBed tool to generate the necessary input BED file. 

**3.2** To process the generated methylation state file (output from Step 3.1) use the methStatesFunction.R script. Relevant functions are explained in comments within the script. The functions allow entropy calculations, epiallele frequency calculations, etc.  An example of running these functions is available in the ExampleUseOfFunctions.R script.


