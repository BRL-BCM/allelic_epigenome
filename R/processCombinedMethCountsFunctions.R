
#this librabry is required for the following script
library(GenomicRanges)

chrOrder = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")


# function that imports the combined methylation counts file into R as GenomicRanges
import_meth_ai_table_combined = function(fileName){
  inputDf = read.table(fileName,header=FALSE,sep="\t",na.strings=c("."))
  order = c();
  print("Done reading!")
  for(chr in chrOrder){
    order = c(order,which(inputDf[,1]==chr))
  }
  inputDf = inputDf[order,]
  print("Done reordering!")
  gRangesInput = GRanges(seqnames = Rle(as.character(inputDf[,1])),
                         ranges = IRanges(start = as.numeric(inputDf[,3]),width = rep(1,nrow(inputDf)),names = paste(inputDf[,1],inputDf[,3],sep=":")),
                         Allele1.Methylated = inputDf[,4],
                         Allele1.Unmethylated = inputDf[,5],
                         Allele2.Methylated = inputDf[,6],
                         Allele2.Unmethylated = inputDf[,7],
                         Total.number.of.reads = inputDf[,8],
                         Maximum.number.of.CpGs = inputDf[,9],
                         Number.of.Heterozygous.Samples = inputDf[,10],
                         Allele1 = as.character(inputDf[,11]),
                         Allele2 = as.character(inputDf[,12]),
                         rsID = inputDf[,13],
                         Allele1.Frequency = inputDf[,14],
                         Allele2.Frequency = inputDf[,15],
                         Ancestral.Allele = inputDf[,16],
                         G1K.Pilot.Mask = inputDf[,17],
                         G1K.Strict.Mask = inputDf[,18],
                         Is.Reference.Allele.Present = inputDf[,19])
  return(gRangesInput)
}

#function for regions of interest
import_ROIs = function(fileName,skip=0,header=TRUE,zeroBased=TRUE){
  inputDf = read.table(fileName,sep="\t",skip=skip,header=header,comment="")
  inputDf = inputDf[-which(!inputDf[,1]%in%chrOrder),]
  gRangesInput=GRanges()
  if(ncol(inputDf)<4){
    inputDf[,4] = seq(1,nrow(inputDf))
    inputDf[,5] = rep(0,nrow(inputDf))
    inputDf[,6] = rep("*",nrow(inputDf))
  }
  if(ncol(inputDf)<5){
    inputDf[,5] = rep(0,nrow(inputDf))
    inputDf[,6] = rep("*",nrow(inputDf))
  }
  if(zeroBased){
    gRangesInput = GRanges(seqnames = Rle(as.character(inputDf[,1])),
                           ranges = IRanges(start = as.numeric(inputDf[,2]+1),end = as.numeric(inputDf[,3]),names = as.character(inputDf[,4])),strand=inputDf[,6])
  } else{
    gRangesInput = GRanges(seqnames = Rle(as.character(inputDf[,1])),
                           ranges = IRanges(start = as.numeric(inputDf[,2]),end = as.numeric(inputDf[,3]),names = as.character(inputDf[,4])),strand=inputDf[,6])    
  }
  gRangesInput = sort(gRangesInput)
  return(gRangesInput)
}


# Compute and add p-values and FDR to genomic ranges using this function
compute_p_values_methylation = function(aiTable,minCounts=6){
  compute_p = function(x){
    p = try(fisher.test(matrix(x,ncol=2))$p.value,silent=TRUE)
    if(is.character(p)){
      warning(p)
      NA
    } else{
      p
    }
  }
  counts = as.matrix(mcols(aiTable)[,1:4])
  P.Values = apply(counts,1,compute_p)
  FDR = p.adjust(P.Values,"fdr")
  mcols(aiTable) = as.data.frame(cbind(as.data.frame(mcols(aiTable)),P.Values,FDR))
  return(aiTable)
}

compute_methylation_al1 = function(aiTable,minCounts=5){
  compute_meth = function(x){
    methLevel = (x[1])/(x[1]+x[2])
    methLevel
  }
  counts = as.matrix(mcols(aiTable)[,1:4])
  goodSnps = which((counts[,1]+counts[,2]>=minCounts))
  Methylation.Allele1 = rep(NA,nrow(counts))
  Methylation.Allele1[goodSnps] = apply(counts[goodSnps,],1,compute_meth)
  mcols(aiTable) = as.data.frame(cbind(as.data.frame(mcols(aiTable)),Methylation.Allele1))
  return(aiTable)
}

compute_methylation_al2 = function(aiTable,minCounts=5){
  compute_meth = function(x){
    methLevel = (x[3])/(x[3]+x[4])
    methLevel
  }
  counts = as.matrix(mcols(aiTable)[,1:4])
  goodSnps = which((counts[,3]+counts[,4]>=minCounts))
  Methylation.Allele2 = rep(NA,nrow(counts))
  Methylation.Allele2[goodSnps] = apply(counts[goodSnps,],1,compute_meth)
  mcols(aiTable) = as.data.frame(cbind(as.data.frame(mcols(aiTable)),Methylation.Allele2))
  return(aiTable)
}


#function to compute allelic methylation differences on genomic ranges object
compute_methylation_differences = function(aiTable,minCounts=5){
  compute_meth_diff = function(x){
    methLevel1 = x[1]/(x[1]+x[2])
    methLevel2 = x[3]/(x[3]+x[4])
    methLevel1-methLevel2
  }
  counts = as.matrix(mcols(aiTable)[,1:4])
  goodSnps = which((counts[,1]+counts[,2]>=minCounts) & (counts[,3]+counts[,4]>=minCounts) & seqnames(aiTable)!="chrY")
  Methylation.Difference = rep(NA,nrow(counts))
  Methylation.Difference[goodSnps] = apply(counts[goodSnps,],1,compute_meth_diff)
  mcols(aiTable) = as.data.frame(cbind(as.data.frame(mcols(aiTable)),Methylation.Difference))
  return(aiTable)
}


#function to export processed combined methylation counts (with the p-values, FDR, and methylation differences added using previous functions)
export_processed_meth_counts = function(methCounts,outFileName){
  methCounts = as.data.frame(methCounts)
  methCounts = methCounts[,-c(4,5)]
  colnames(methCounts)[1] = "chromosome"
  methCounts[,2] = methCounts[,2]-1
  write.table(methCounts,file = outFileName,quote=FALSE,sep="\t",row.names = FALSE)
}
