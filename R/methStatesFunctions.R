aiDistPerStability <- function(aisPerRegion,stability){
  maxNumStates <- max(stability)
  result <- NULL
  for(i in 1:maxNumStates){
    regsPerStability <- unique(gsub("\\.window.+$","",names(stability)[stability==i]))
    result <- rbind(result,
                    c(sum(aisPerRegion[regsPerStability]==0,na.rm=T),
                    sum(aisPerRegion[regsPerStability]>=1,na.rm=T)))
  }
  result <- cbind(result,result[,2]/rowSums(result))
  rownames(result) <- paste("nStates",1:maxNumStates,sep="=")
  colnames(result) <- c("noAI","withAI","propWithAI")
  return(result)
}

findStability <- function(patternMatrix,propOfObs=0.05){
  findNumStates <- function(stateDist){
    #min(which(cumsum(sort(stateDist/sum(stateDist),decreasing=T))>=propOfObs)))
    return(sum((stateDist/sum(stateDist))>=propOfObs))
  }
  result <- apply(patternMatrix,1,findNumStates)
  return(result)
}


# function takes pattern and entropies output from function extractPatternsandEntropiesAllRegions(), window size (should be same as specified for extractPatternsandEntropiesAllRegions function), minimum required Coverage (e.g. 30 used in past), minimum entropy, and maximum entropy as arguments. Output is matrix with all possible states as columns and rows are windows in that region. If a region or window didn’t meet criteria is gets dropped from output matrix.

getObservedPatternsMatrix <- function(patAndEnt,windowSize,minCoverage,minEntropy,maxEntropy){
  resultVec <- rep(0,2^windowSize)
  names(resultVec) <- getPossibleStates(windowSize)
  result <- lapply(patAndEnt,function(region){
              lapply(region,function(window){
                if(is.null(window)){return(NULL)}
                if(window$total.num.obs>=minCoverage & 
                   window$entropy <= maxEntropy &
                   window$entropy >= minEntropy){
                  resultVec[colnames(window$obs.per.state)] <- window$obs.per.state
                  return(resultVec)
                } else{
                  return(NULL)
                }
              })
            })
  result <- unlist(result)
  windowNames <- unique(gsub("(region\\.\\d+\\.window\\.\\d+)(\\.\\d+$)","\\1",names(result),perl=T))
  result <- matrix(result,ncol=2^windowSize,byrow=T)
  colnames(result) <- names(resultVec)
  rownames(result) <- windowNames
  return(result)
}

computeMeanSlidingEntropyPerRegion <- function(methStates,windowSize=3){
  computeEntropyOneRegion <- function(methStatesForRegion){
    if(is.null(methStatesForRegion$patterns)){
      return(c(NA,NA))
    } else{
      return(computeMeanSlidingEntropy(methStatesForRegion$patterns,windowSize))
    }
  }
  result <- matrix(unlist(lapply(methStates,computeEntropyOneRegion)),ncol=2,byrow=T)
  colnames(result) <- c("Entropy","Max.Entropy")
  rownames(result) <- names(methStates)
  return(result)
}


# function takes the variable generated from previous function (readMethStates()),window size (number of CpGs user want to group to test e.g. for 4 CpG window = 4), and minimum coverage as parameters. Output is list that has one element per region. Each element has each one of the sliding windows. Inside each window is information for that window (entropy, max entropy, total # observations (total number of reads that mapped there and had information for all CpGs in that window (if at least one 0 value for a CpG in window in that read, it doesn’t count), and observations per state). If a region didn’t have any window meeting the specified criteria, will be returned NULL.

extractPatternsAndEntropiesAllRegions <- function(methStates,windowSize=3,minCoverage=0){
  computeEntropyOneRegion <- function(methStatesForRegion){
    if(is.null(methStatesForRegion$patterns)){return(NULL)}
    computeAllSlidingEntropies(methStatesForRegion$patterns,windowSize,minCoverage)
  }
  result <- lapply(methStates,computeEntropyOneRegion)
  return(result)
}

computeAllSlidingEntropies <- function(patterns,windowSize,minCoverage=0){
  fullPatternSize <- nchar(patterns[1,1])
  if(fullPatternSize < windowSize){
    return(NULL)
  }
  computeEntropyForWindow <- function(window){
    temp <- computeEntropy(substr(patterns[,1],window[1],window[2]),patterns[,3],minCoverage)
    if(is.null(temp)){return(NULL)}
    return(list(entropy=temp$entropy,
                max.entropy=temp$max.entropy,
                total.num.obs=temp$total.num.obs,
                obs.per.state=temp$obs.per.state))
  }
  windows <- 1:(fullPatternSize-windowSize+1)
  windows <- cbind(windows,windows+(windowSize-1))
  rownames(windows) <- paste("window",1:nrow(windows),sep=".")
  result <- apply(windows,1,computeEntropyForWindow)
  return(result)
}

computeMeanSlidingEntropy <- function(patterns,windowSize){
  fullPatternSize <- nchar(patterns[1,1])
  if(fullPatternSize < windowSize){
    return(c(NA,NA))
  }
  computeEntropyForWindow <- function(window){
    temp <- computeEntropy(substr(patterns[,1],window[1],window[2]),patterns[,3])
    if(is.null(temp)){return(c(NA,NA))}
    return(c(temp$entropy,temp$max.entropy))
  }
  windows <- 1:(fullPatternSize-windowSize+1)
  windows <- cbind(windows,windows+(windowSize-1))
  entropies <- apply(windows,1,computeEntropyForWindow)
  rownames(entropies) <- c("Entropy","Max.Entropy")
  if(fullPatternSize>=windowSize+1){
    entropies <- rowMeans(entropies,na.rm=T)
  } else{
    entropies <- entropies[,1]
  }
  return(entropies)
}

computeEntropy <- function(cpgPatterns,counts,minCoverage=0){
  expandPatterns <- function(x){
    return(rep(cpgPatterns[x],counts[x]))
  }
  cpgPatterns <- unlist(sapply(1:length(cpgPatterns),expandPatterns))
  cpgPatterns <- sort(grep("0",cpgPatterns,value=T,invert=T))
  if(length(cpgPatterns)<=minCoverage){return(NULL)}
  
  obsPerState <- t(as.matrix(table(cpgPatterns)))
  numObs <- length(cpgPatterns)
  entropy <- -sum((obsPerState/numObs) * log2(obsPerState/numObs))
  numPossibleStates <- 2^nchar(cpgPatterns[1])
  
  maxEntropy <- -numPossibleStates*((1/numPossibleStates)*log2(1/numPossibleStates))
  if(numObs < numPossibleStates){
    maxEntropy <- -numObs*((1/numObs)*log2(1/numObs))
  }
  
  return(list(entropy = entropy,
              max.entropy = maxEntropy,
              num.possible.states = numPossibleStates,
              num.obs.states = length(obsPerState),
              total.num.obs = numObs,
              obs.per.state = obsPerState))
}


# reads in meth state file generated in step 3.1 (see ReadMe file). Returns a list of regions of interest.Each region of interest has a list that has the columns: chr,start, patterns, etc. (all info is for the + strand) Patterns are a matrix that has 3 columns: 1st column has pattern of the CpGs (0 = no info, 1 = Unmethylated, 2 = methylated), 2nd column has info about the hetSNP (* = not available), 3rd is counts for each observations of the patterns of het:CpG pattern pair).
readMethStates <- function(file){
  inputDf <- read.table(file,header=T)
  processInput <- function(x){
    chr <- as.character(x[1])
    start <- as.numeric(as.character(x[2]))
    end <- as.numeric(as.character(x[3]))
    cpg.pos <- as.numeric(unlist(strsplit(as.character(x[4]),split = ",")))
    var.pos = NULL
    if(x[5] != "."){
      var.pos <- as.numeric(unlist(strsplit(as.character(x[5]),split = ",")))
    }
    patterns = NULL
    if(x[6] != "."){
      temp <- unlist(strsplit(as.character(x[6]),split = ","))
      patterns <- matrix(unlist(strsplit(temp,":")),nrow=length(temp),ncol=3,byrow=T)
      patterns <- as.data.frame(patterns)
      patterns[,1] <- as.character(patterns[,1])
      if(is.null(var.pos)){
        patterns[,2] = rep(NA,nrow(patterns))
      } else{
        patterns[,2] <- as.character(patterns[,2])
      }
      patterns[,3] <- as.numeric(as.character(patterns[,3]))
    }
    return(list(chr = chr,
                start = start,
                end = end,
                cpg.pos = cpg.pos,
                var.pos = var.pos,
                patterns = patterns))
  }
  out <- apply(inputDf,1,processInput)
  names(out) <- paste("region",1:nrow(inputDf),sep=".")
  return(out)
}

getPossibleStates <- function(windowSize){
  possibleStates <- c("1","2")
  temp <- NULL
  for(i in 2:windowSize){
    for(s in possibleStates){
      s1 <- paste(s,"1",sep="")
      s2 <- paste(s,"2",sep="")
      temp <- c(temp,s1,s2)
    }
    possibleStates <- temp
    temp <- NULL
  }
  return(possibleStates)
}
