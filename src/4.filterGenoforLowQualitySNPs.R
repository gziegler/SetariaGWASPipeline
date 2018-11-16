library(data.table)
load("../data/genotype/3.FilteredGenotypeFile.MatrixFormat.noscaffold.hetFilter0.25.maf0.1.rda")

dist <- 1
neighborCors <- sapply(1:(nrow(genoMatrix)-dist),function(x) {if(x%%100000==0){print(x)};cor(genoMatrix[x,],genoMatrix[x+dist,],use="complete.obs")})
neighborCors <- neighborCors^2

#neighborDists <- sapply(1:(nrow(alleleTable)-dist),function(x){alleleTable$bp[x+dist]-alleleTable$bp[x]})
neighborDists <- alleleTable$pos[2:nrow(alleleTable)]-alleleTable$pos[1:(nrow(alleleTable)-1)]
#now we get the correct number (8) of less than 0 distances indicating break between chroms, set those to NA
neighborDists[neighborDists<0] <- NA
summary(neighborDists)


pdf("../results/4.FilteredGeno.neighborCor.pdf")
hist(neighborCors,breaks=100)
dev.off()

neighbors <- data.table(chr=alleleTable$chr[-nrow(alleleTable)],pos=alleleTable$pos[-nrow(alleleTable)],neighborDists=neighborDists,neighborCors=neighborCors)


pdf("../results/4.FilteredGeno.ChromsomeWideNeighboringSNP.LD.subset.pdf",width=20)
for(i in 1:9){
  subsetN <- neighbors[neighbors$chr==i]
  subsetN <- subsetN[sample(1:nrow(subsetN),size = 10000,replace = FALSE)]
  plot(neighborCors~pos,data=subsetN,type="p", ylim=c(0,1),pch=20,xlim=c(1,max(pos)),col="navy",ylab="r^2",cex.lab=1.6,main=paste("Chromosome",i))
}
dev.off()


#remove non-correlated SNPs
###Find a good SNP call (e.g. from end of stretch of highly correlated)
####look at SNP and it's neighbor, if correlated, move to next, if not, delete next, and move to next after that
#####check stop criteria as well
#return indices of genotype/snpInfo to keep
#neighborCors=list of precomputed correlations
#highCorThreshold=threshold for finding 'good' SNPs defined as groupings of 2 or more SNPs where correlation between each neighbor is > threshold
#breakDist=base pair distance to break to the next 'good' SNP
#deleteThreshold=delete SNP if not correlated at > threshold with current 'good' snp unless > breakDist away
filterSNPsbyCor <- function(genotype,snpInfo,neighborCors,highCorThreshold=0.9,deleteThreshold=0.5,breakDist=2000){
  dist <- 1 #distance for calculating neighborCors
  if(missing(neighborCors)){
    cat("Calculating Neighbor Correlations")
    neighborCors <- sapply(1:(nrow(genotype)-dist),function(x) {if(x%%100000==0){cat(".")};cor(genotype[x,],genotype[x+dist,],use="complete.obs")})
    neighborCors <- neighborCors^2
    print("Done")
  }
  
  #   if(missing(neighborDists)){
  #     cat("Calculating Neighbor Distances")
  #     neighborDists <- sapply(1:(nrow(snpInfo)-dist),function(x){snpInfo$bp[x+dist]-snpInfo$bp[x]})
  #     print("Done")
  #   }
  
  highCorIdxs <- which(neighborCors>=highCorThreshold)
  s <- split(highCorIdxs, cumsum(c(0, diff(highCorIdxs) != 1)))
  highCorClusters <- unname(s[which(unname(sapply(s, length))>1)])
  genoIndx <- unlist(highCorClusters[1])[1]
  highCorIdx <- 2 #index of next 'good' SNP cluster
  backSearchIndx <- 1 #index of last SNP tested when skipped to next 'good' SNP
  forwardSearchIndx <- genoIndx+1 #index of next SNP to be tested during forward search
  keepIdxs <- genoIndx
  possibleMissedSNPs <- c()
  cat("Starting removal of bad SNPs:")
  while((genoIndx<nrow(genotype) & forwardSearchIndx<nrow(genotype))){
    if(length(keepIdxs)%%20000==0){print(paste(length(keepIdxs),"SNPs kept at SNP number",genoIndx,"out of",nrow(genotype),"SNPs.",length(possibleMissedSNPs),"possible missed SNPs."))}
    #backward search
    if(backSearchIndx<genoIndx){
      preBackIndx <- genoIndx
      for(i in (genoIndx-1):backSearchIndx){
        thisCor <- ifelse((genoIndx-i)==1,neighborCors[i],cor(genotype[i,],genotype[genoIndx,],use="complete.obs")^2)
        if(thisCor >= deleteThreshold){
          genoIndx <- i
          keepIdxs <- c(keepIdxs,i)
        }else if((snpInfo$pos[genoIndx]-snpInfo$pos[i])>=breakDist){
          #print(paste("Could not find 'good' SNP within",breakDist,"bp threshold during backward search. Comparing it to",snpInfo$snpNames[genoIndx],"Deleting SNP",snpInfo$snpNames[i]))
          possibleMissedSNPs <- c(possibleMissedSNPs,paste(snpInfo$chr[i],snpInfo$pos[i],sep="_"))
        }
      }
      backSearchIndx <- preBackIndx
      genoIndx <- preBackIndx
      rm(thisCor,preBackIndx)
    }
    
    #forward search
    thisCor <- ifelse((forwardSearchIndx-genoIndx)==1,neighborCors[genoIndx],cor(genotype[genoIndx,],genotype[forwardSearchIndx,],use="complete.obs")^2)
    if(thisCor >= deleteThreshold){
      genoIndx <- forwardSearchIndx
      backSearchIndx <- genoIndx
      forwardSearchIndx <- genoIndx+1
      keepIdxs <- c(keepIdxs,genoIndx)
    }else if((snpInfo$pos[forwardSearchIndx]-snpInfo$pos[genoIndx])>=breakDist){ #snp doesn't meet threshold, is it outside bp distance
      prevBackSearchIndx <- backSearchIndx
      backSearchIndx <- forwardSearchIndx
      #######test to make sure didn't move through a high cor cluster#####
      if(backSearchIndx > unlist(highCorClusters[highCorIdx])[1]){
        while(backSearchIndx > unlist(highCorClusters[highCorIdx])[1]){
          highCorIdx <- highCorIdx + 1
          if(highCorIdx>length(highCorClusters)){
            highCorIdx <- length(highCorClusters)
            print("Caution: Couldn't find 'good' SNP at end of genotype, increasing break distance to sample to end of chromosome.")
            breakDist <- 1e6
            break
          }
        }
      }
      if(unlist(highCorClusters[highCorIdx])[1]>genoIndx){
        genoIndx <- unlist(highCorClusters[highCorIdx])[1]
        keepIdxs <- c(keepIdxs,genoIndx)
      }
      forwardSearchIndx <- genoIndx+1
    }else{
      forwardSearchIndx <- forwardSearchIndx+1
    }
  }#end while
  list(keeps=sort(keepIdxs),possibleMissed=possibleMissedSNPs)
}



filterResults <- list()
filterGeno <- matrix()
filterInfo <- data.table()
for(i in 1:9){
  print(paste("Running on chromsome",i))
  thisNeighborCor <- neighborCors[which(alleleTable$chr==i)]
  thisGeno <- genoMatrix[which(alleleTable$chr==i),]
  thisSnpInfo <- alleleTable[which(alleleTable$chr==i),]
  filterOut <- filterSNPsbyCor(thisGeno,thisSnpInfo,thisNeighborCor)
  if(nrow(filterGeno)<=1){
    filterGeno <- thisGeno[filterOut$keeps,]
  }else{
    filterGeno <- rbind(filterGeno,thisGeno[filterOut$keeps,])
  }
  filterInfo <- rbind(filterInfo,thisSnpInfo[filterOut$keeps,])
  filterResults[[i]] <- filterOut
  rm(thisNeighborCor,thisGeno,thisSnpInfo,filterOut)
}

sum(sapply(filterResults,function(x)length(x$keeps)))

dist <- 1
filterCors <- sapply(1:(nrow(filterGeno)-dist),function(x) {if(x%%100000==0){cat(".")};cor(filterGeno[x,],filterGeno[x+dist,],use="complete.obs")})
filterCors <- filterCors^2

filterDists <- diff(filterInfo$pos)
filterDists[filterDists<0] <- NA
summary(filterDists)

neighbors <- data.table(chr=filterInfo$chr[-nrow(filterInfo)],pos=filterInfo$pos[-nrow(filterInfo)],neighborDists=filterDists,neighborCors=filterCors)


pdf("../results/4.BadSNPsFiltered.ChromsomeWideNeighboringSNP.LD.subset.pdf",width=20)
for(i in 1:9){
  subsetN <- neighbors[neighbors$chr==i]
  subsetN <- subsetN[sample(1:nrow(subsetN),size = 10000,replace = FALSE)]
  plot(neighborCors~pos,data=subsetN,type="p", ylim=c(0,1),pch=20,xlim=c(1,max(pos)),col="navy",ylab="r^2",cex.lab=1.6,main=paste("Chromosome",i))
}
dev.off()

save(filterInfo,filterGeno,neighbors,file="../data/genotype/4.filteredSNPs.2kbDistThresh.0.5neighborLD.rda")

