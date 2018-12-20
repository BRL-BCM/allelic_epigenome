# init represesnts a .txt file (init.txt) created by the user that contains the following 12 columns (include these headers in your file): chromosome, start, end, Allele1.Methylated, Allele1.Unmethylated, Allele2.Methylated, Allele2.Unmethylated, Total.number.of.reads, Maximum.number.of.CpGs, Number.of.Heterozygous.Samples, Allele1, Allele2
# Columns 1-3 should be filled with all of the  loci of interest where allelic imbalance calling for methylation will be performed (e.g. chr1 529781 529782)
# Columns 4-10 should be filled with values of 0 as placeholders for all of these counts
# Columns 11-12 should be filled with the two alleles present at the locus (e.g A T)
# Note: Edit the path and file name (if the init.txt file was named differently) to the directory where the init.txt file is stored


init = read.table("/mnt/brlstor/cluster.data2/allelicEpigenome/methCounts_v2/combined/init.txt.gz",header = T, sep = '\t',stringsAsFactors = FALSE)


# samples contains the list of methylation counts files to be combined. Change the path to the directory where the individual allelic methylation count files (output from Step 2.2) are located

samples = list.files(path = "/cluster.data2/allelicEpigenome/methCounts_v2/autosome/",full.names = T,recursive = FALSE)


# A for loop that loops through the individual methylation counts files present in the samples variable to commbine counts for both alleles (both methylated and unmethylated), combine the total number of reads used for the counts across the samples, the maximum number of CpGs surrounding the locus, and the total number of individual samples that were in a heterozygous state

for(i in samples){
  new = read.table(i,header = F, sep = '\t',skip =1,stringsAsFactors = FALSE)
  colnames(new) = c('chromosome','start','end',"A1_M","A1_U","A2_M","A2_U","reads","CpGs","A1","A2")
  init = merge(init,new,by = c('chromosome',"start","end"),all.x = T)
  init[is.na(init)] <- 0
  init$Total.number.of.reads = init$Total.number.of.reads+init$reads
  init$Maximum.number.of.CpGs = ifelse(init$Maximum.number.of.CpGs >= init$CpGs,init$Maximum.number.of.CpGs,init$CpGs)
  init$Number.of.Heterozygous.Samples = ifelse(init$A1 == init$A2,init$Number.of.Heterozygous.Samples,init$Number.of.Heterozygous.Samples+1)
  init$Allele1.Methylated = with(init,ifelse(init$Allele1 == init$A1 & init$Allele1 == init$A2,init$Allele1.Methylated+init$A1_M+init$A2_M,ifelse(init$Allele1 == init$A1,init$Allele1.Methylated+init$A1_M,ifelse(init$Allele1 == init$A2,init$Allele1.Methylated+init$A2_M,init$Allele1.Methylated))))
  init$Allele1.Unmethylated = with(init,ifelse(init$Allele1 == init$A1 & init$Allele1 == init$A2,init$Allele1.Unmethylated+init$A1_U+init$A2_U,ifelse(init$Allele1 == init$A1,init$Allele1.Unmethylated+init$A1_U,ifelse(init$Allele1 == init$A2,init$Allele1.Unmethylated+init$A2_U,init$Allele1.Unmethylated))))
  init$Allele2.Methylated = with(init,ifelse(init$Allele2 == init$A1 & init$Allele2 == init$A2,init$Allele2.Methylated+init$A1_M+init$A2_M,ifelse(init$Allele2 == init$A1,init$Allele2.Methylated+init$A1_M,ifelse(init$Allele2 == init$A2,init$Allele2.Methylated+init$A2_M,init$Allele2.Methylated))))
  init$Allele2.Unmethylated = with(init,ifelse(init$Allele2 == init$A1 & init$Allele1 == init$A2,init$Allele2.Unmethylated+init$A1_U+init$A2_U,ifelse(init$Allele2 == init$A1,init$Allele2.Unmethylated+init$A1_U,ifelse(init$Allele2 == init$A2,init$Allele2.Unmethylated+init$A2_U,init$Allele2.Unmethylated))))
  init = init[1:12]
  rm(new)
  }


# edit path and file name to output the combined methCount file into the desired directory

write.table(init,"/cluster.data2/allelicEpigenome/methCounts_v2/autosome/combinedMethCounts_v2.txt",col.names = T,row.names = F,sep ='\t',quote = F)


