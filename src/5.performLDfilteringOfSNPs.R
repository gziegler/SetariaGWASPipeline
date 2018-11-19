library(data.table)
load("../data/genotype/4.filteredSNPs.2kbDistThresh.0.5neighborLD.rda")

#Filter SNPs with correlation >0.975
#filter SNPs
#recalc neighbor cors
#if still neighbor cors >0.975, repeat
filterSNPsinHighCor <- function(genotype,neighborCors,threshold,dist=1){
  #  done <- FALSE
  #  while(!done){
  if(missing(neighborCors)){
    cat("Calculating Neighbor Correlations")
    neighborCors <- sapply(1:(nrow(genotype)-dist),function(x) {if(x%%100000==0){cat(".")};cor(genotype[x,],genotype[x+dist,],use="complete.obs")})
    neighborCors <- neighborCors^2
    print("Done")
  }
  #testIdxs <- which(neighborCors>threshold)
  #s <- split(testIdxs, cumsum(c(0, diff(testIdxs) != 1)))
  #highCorClusters <- unname(s[which(unname(sapply(s, length))>1)])
  #print("Filtering",length(testIdxs),"SNP Pairs")
  #    sapply(testIdxs,functions(x){x})
  genoIndx <- 1
  forwardSearchIndx <- genoIndx+1
  keepIdxs <- genoIndx
  cat(paste0("Starting removal of highly correlated SNPs with correlation >",threshold,":"))
  while((genoIndx<nrow(genotype) & forwardSearchIndx<nrow(genotype))){
    if(length(keepIdxs)%%20000==0){print(paste(length(keepIdxs),"SNPs kept at SNP number",genoIndx,"out of",nrow(genotype),"SNPs."))}
    thisCor <- ifelse((forwardSearchIndx-genoIndx)==1,neighborCors[genoIndx],cor(genotype[genoIndx,],genotype[forwardSearchIndx,],use="complete.obs")^2)
    #if SNP is in high LD with neighbor, keep it, if not, move genoIndx
    if(thisCor >= threshold){
      forwardSearchIndx <- forwardSearchIndx+1
    }else{
      genoIndx <- forwardSearchIndx
      forwardSeachIndx <- genoIndx+1
      keepIdxs <- c(keepIdxs,genoIndx)
    }
  }#end while through genotype
  keepIdxs
  #update genotype, recalc cors, decide if we need to do another pass
  #genotype <- genotype[keepIdxs,]
  #snpInfo <- snpInfo[keepIdxs,]
  #}#end while through iterations
}



filterHighResults <- list()
filterHighGeno <- matrix()
filterHighInfo <- data.table()
for(i in 1:9){
  print(paste("Running on chromsome",i))
  thisNeighborCor <- neighbors$neighborCors[which(neighbors$chr==i)]
  thisGeno <- filterGeno[which(filterInfo$chr==i),]
  thisSnpInfo <- filterInfo[which(filterInfo$chr==i),]
  filterOut <- filterSNPsinHighCor(thisGeno,thisNeighborCor,0.975)
  if(nrow(filterHighGeno)<=1){
    filterHighGeno <- thisGeno[filterOut,]
  }else{
    filterHighGeno <- rbind(filterHighGeno,thisGeno[filterOut,])
  }
  filterHighInfo <- rbind(filterHighInfo,thisSnpInfo[filterOut,])
  filterHighResults[[i]] <- filterOut
  rm(thisNeighborCor,thisGeno,thisSnpInfo,filterOut)
}



dist <- 1
filterHighCors <- sapply(1:(nrow(filterHighGeno)-dist),function(x) {if(x%%100000==0){cat(".")};cor(filterHighGeno[x,],filterHighGeno[x+dist,],use="complete.obs")})
filterHighCors <- filterHighCors^2

filterHighDists <- diff(filterHighInfo$pos)
filterHighDists[filterHighDists<0] <- NA
summary(filterHighDists)

neighbors <- data.table(chr=filterHighInfo$chr[-nrow(filterHighInfo)],pos=filterHighInfo$pos[-nrow(filterHighInfo)],neighborDists=filterHighDists,neighborCors=filterHighCors)

pdf("../results/5.FilteredHighSNP.ChromsomeWideNeighboringSNP.LD.subset.pdf",width=20)
for(i in 1:9){
  subsetN <- neighbors[neighbors$chr==i]
  subsetN <- subsetN[sample(1:nrow(subsetN),size = 10000,replace = FALSE)]
  plot(neighborCors~pos,data=subsetN,type="p", ylim=c(0,1),pch=20,xlim=c(1,max(pos)),col="navy",ylab="r^2",cex.lab=1.6,main=paste("Chromosome",i))
}
dev.off()

#xz compression is super slow to save (<10min), but compresses enough (81mb vs 120mb) to get file size below 100mb for github and still uncompresses about as fast (30s)
save(neighbors,filterHighInfo,filterHighGeno,file="../data/genotype/5.filteredSNPs.noHighCorSNPs.2kbDistThresh.0.5neighborLD.0.975LDfilter.rda",compress="xz",compression_level = 5)

