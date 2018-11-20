###1) further thin genotype file by LD to get to better number for structure/kinship calculation
#######The median value is 0.78. sp isomg 0.7 for cutoff will get rid of about 60% of the SNPs
###2) Calculate PCs for pop structure correction
###3) Calculate Kinship matrix using Astle-Balding method
rm(list=ls())
if(!require(ionomicsUtils,quietly = TRUE)){
  library(devtools)
  install_github("gziegler/ionomicsUtils")
}
library(ionomicsUtils)
library(data.table)
library(gplots)

load("../data/genotype/5.filteredSNPs.noHighCorSNPs.2kbDistThresh.0.5neighborLD.0.975LDfilter.rda")

#Filter SNPs with correlation >theshold
#filter SNPs
#recalc neighbor cors
#if still neighbor cors >threshold, repeat
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

structFiltResults <- list()
structFiltGeno <- matrix()
structFiltInfo <- data.table()
for(i in 1:9){
  print(paste("Running on chromsome",i))
  thisNeighborCor <- neighbors$neighborCors[which(neighbors$chr==i)]
  thisGeno <- filterHighGeno[which(filterHighInfo$chr==i),]
  thisSnpInfo <- filterHighInfo[which(filterHighInfo$chr==i),]
  filterOut <- filterSNPsinHighCor(thisGeno,thisNeighborCor,0.7)
  if(nrow(structFiltGeno)<=1){
    structFiltGeno <- thisGeno[filterOut,]
  }else{
    structFiltGeno <- rbind(structFiltGeno,thisGeno[filterOut,])
  }
  structFiltInfo <- rbind(structFiltInfo,thisSnpInfo[filterOut,])
  structFiltResults[[i]] <- filterOut
  rm(thisNeighborCor,thisGeno,thisSnpInfo,filterOut)
}

if(nrow(structFiltGeno)>250000){
  set.seed(129387)
  keepIdxs <- sort(sample(1:nrow(structFiltGeno),250000,replace = FALSE))
  structFiltGeno <- structFiltGeno[keepIdxs,]
  structFiltInfo <- structFiltInfo[keepIdxs,]
}

genoMat <- t(structFiltGeno)
colnames(genoMat) <- paste(structFiltInfo$chr,structFiltInfo$pos,sep="_")
#genoMat[which(is.na(genoMat))] <- 0

kinship <- ionomicsUtils::realizedAB(genoMat,maf=structFiltInfo$MAF)
dir.create(file.path(paste0("../results")), showWarnings = FALSE)
pdf(paste("../results/6.AstleBalding.synbreed.kinship.pdf",sep=""), width = 12, height = 12)
par(mar = c(25,25,25,25))
heatmap.2(kinship,  cexRow =.2, cexCol = 0.2, col=rev(heat.colors(256)), scale="none", symkey=FALSE, trace="none")
dev.off()

save(kinship,file="../data/genotype/6.AstleBalding.synbreed.kinship.rda")


####Genotype file should be in 0, 1, 2 coding
###rows are snps and columns are individuals
e<-ionomicsUtils::eigenstrat(structFiltGeno)
pdf("../results/6.Eigenstrat.structure.PC1vsPC2.and.VarianceExplained.pdf")
plot(e$vectors[,1:2],col=1,xlab="PC1",ylab="PC2")
barplot(height=e$values[1:20]/sum(e$values))
dev.off()
structData <- e$vectors[,1:50]
rownames(structData) <- colnames(structFiltGeno)
save(structData,file="../data/genotype/6.Eigenstrat.population.structure.50PCs.rda")
write.table(structData,file="../data/genotype/6.Eigenstrat.population.structure.50PCs.csv",sep=",",col.names=TRUE,row.names=TRUE)
write.table(e$values[1:50]/sum(e$values),file="../results/6.Eigenstrat.population.structure.varianceExplained.csv",sep=",",col.names=TRUE,row.names=TRUE)
