library(feather)
library(data.table)
library(plyr)
library(reshape2)
library(dplyr)
resultsFiles <- list.files("../GWASresults/clumpRes","\\.feather",full.names=TRUE)
dataList <- list()
for(i in resultsFiles){
  message("Reading results file ",i,"\n")
  dataList[[length(dataList)+1]] <- feather::read_feather(path=i)
}
snpTable <- data.table::rbindlist(dataList)
rm(dataList)


pvalCols <- grep("pval|fisherP",colnames(snpTable),value=TRUE)

summariseSNP <- snpTable[, lapply(.SD, min,na.rm=TRUE),by=clump,.SDcols = pvalCols] #get minimum pvalue for each snp from each clump

#Add chromsome and range to each SNP
clumpPos <- snpTable %>% group_by(clump) %>% dplyr::summarise(chrom=min(chrom),minBP=min(pos),maxBP=max(pos),nSNPs=dplyr::n())
summariseSNP <- merge(clumpPos,summariseSNP,by="clump")
setDT(summariseSNP)
#recalculate fisher pvals
#Deconvulte Trait ID into a table
pvalCols <- grep("pval",colnames(snpTable),value=TRUE)
traitTable <- data.frame(Phenotype = sapply(strsplit(pvalCols,"_\\s*(?=[^_]+$)",perl=TRUE),"[",1), #regex to split on last _ (sometimes phenotypes have an '_')
                         DAP =       as.numeric(sapply(strsplit(sapply(strsplit(sapply(strsplit(pvalCols,"_\\s*(?=[^_]+$)",perl=TRUE),"[",2),"\\."),"[",3),"-"),"[",1)),
                         Treatment = sapply(strsplit(sapply(strsplit(pvalCols,"_\\s*(?=[^_]+$)",perl=TRUE),"[",2),"\\."),"[",2),
                         Experiment = sapply(strsplit(sapply(strsplit(pvalCols,"_\\s*(?=[^_]+$)",perl=TRUE),"[",2),"\\."),"[",1),
                         colID = pvalCols,
                         stringsAsFactors = FALSE)

for(pheno in unique(traitTable$Phenotype)){
  for(dap in unique(traitTable$DAP[which(traitTable$Phenotype==pheno)])){
    for(treat in unique(traitTable$Treatment[which(traitTable$Phenotype==pheno & traitTable$DAP==dap)])){
      message(pheno," ",dap," ",treat,"\n")
      thisTraitCols <- grep("pval",traitTable$colID[which(traitTable$Phenotype==pheno & traitTable$DAP==dap & traitTable$Treatment==treat)],value=TRUE)
      chiSqDF = 2*length(thisTraitCols)
      #thisTraitCols <- grep("mean",thisTraitCols,value=TRUE,invert = TRUE) #exclude columns from rank.mean experiment
      summariseSNP[,paste(pheno,treat,dap,"fisherPclump",sep="_") := pchisq((-2*rowSums(log(.SD),na.rm=T)),df=chiSqDF,lower.tail=F), .SDcols = thisTraitCols]
      #mergeRes[,paste(pheno,treat,dap,"fisherP",sep="_") := fisher_combine_pvals(.SD),by=seq_len(nrow(mergeRes)),.SDcols = thisTraitCols]
    }
  }
}

clumpCols <- grep("fisherPclump",colnames(summariseSNP),value=TRUE)
finalThresh <- 1e-5
filterRes <- summariseSNP[summariseSNP[, Reduce(`|`, lapply(.SD, `<=`, finalThresh)),.SDcols = clumpCols]] #keep phenotype columns where fisherP <= finalThresh

###Use this when there are more than one trait
summariseMelt <- melt(filterRes[,c("clump","chrom","minBP","maxBP",clumpCols),with=FALSE],id.vars = c("clump","chrom","minBP","maxBP"))
summariseMelt$fisherLogP <- -log10(summariseMelt$value)
filteredSNP <- summariseMelt[summariseMelt$fisherLogP >= 7,]
filteredSNP$Phenotype = sapply(strsplit(sapply(strsplit(sapply(strsplit(as.character(filteredSNP$variable),"_\\s*(?=[^_]+$)",perl=TRUE),"[",1),"_\\s*(?=[^_]+$)",perl=TRUE),"[",1),"_\\s*(?=[^_]+$)",perl=TRUE),"[",1) #regex to split on last _ (sometimes phenotypes have an '_')
filteredSNP$DAP = as.numeric(sapply(strsplit(sapply(strsplit(as.character(filteredSNP$variable),"_\\s*(?=[^_]+$)",perl=TRUE),"[",1),"_\\s*(?=[^_]+$)",perl=TRUE),"[",2))
filteredSNP$Treatment = sapply(strsplit(sapply(strsplit(sapply(strsplit(as.character(filteredSNP$variable),"_\\s*(?=[^_]+$)",perl=TRUE),"[",1),"_\\s*(?=[^_]+$)",perl=TRUE),"[",1),"_\\s*(?=[^_]+$)",perl=TRUE),"[",2)
filteredSNP$clumpCenter = round((filteredSNP$minBP+filteredSNP$maxBP)/2)
filteredSNP$Experiment = "IB005toIB008"

write.table(filterRes,file="../GWASresults/12.clumpColumn.GWASresults.fisherCombinedPval.LTtraits.IB005toIB008.csv",sep=",",col.names=TRUE,row.names=FALSE)
write.table(filteredSNP,"../results/12.clumpColumn.GWASresults.filtered.ZbrowseReady.IB005toIB008.csv",sep=",",col.names=TRUE,row.names=FALSE)
