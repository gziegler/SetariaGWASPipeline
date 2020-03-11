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
  #add max of previous clump size to current clump
  if(length(dataList) != 1){
    dataList[[length(dataList)]]$clump <- dataList[[length(dataList)]]$clump+max(dataList[[length(dataList)-1]]$clump)
  }
}
snpTable <- data.table::rbindlist(dataList)
rm(dataList)

#Open and merge IB0005 with snpTable clump information
IB005 <- feather::read_feather(path="../GWASresults/11.allLTphenos.IB005to8.CombinedMLMMResultsWithFisherPval.allSNPs.feather")
IB005merge <- merge(IB005,snpTable[,c("SNP","clump","chrom","pos")],by="SNP",all=T)
# #Add chromsome and range to each SNP
# IB005merge$chrom <- as.numeric(sapply(strsplit(IB005merge$SNP,"_"),"[",1))
# IB005merge$pos <- as.numeric(sapply(strsplit(IB005merge$SNP,"_"),"[",2))
rm(IB005)
feather::write_feather(IB005merge,path="../GWASresults/13.clumpColumn.IB005to8.CombinedMLMMResultsWithFisherPval.allSNPs.feather")

pvalCols <- grep("pval|fisherP",colnames(snpTable),value=TRUE)
setDT(snpTable)
summariseSNP <- snpTable[, lapply(.SD, min,na.rm=TRUE),by=clump,.SDcols = pvalCols] #get minimum pvalue for each snp from each clump



clumpPos <- snpTable %>% group_by(clump) %>% dplyr::summarise(chrom=min(chrom),minBP=min(pos),maxBP=max(pos),nSNPs=dplyr::n())
summariseSNP <- merge(clumpPos,summariseSNP,by="clump")
setDT(summariseSNP)
#recalculate fisher pvals
#Deconvulte Trait ID into a table
pvalCols <- grep("pval",colnames(snpTable),value=TRUE)
traitTable <- data.frame(Phenotype = sapply(strsplit(sapply(strsplit(pvalCols,"_\\s*(?=[^_]+$)",perl=TRUE),"[",1),"_\\s*(?=[^_]+$)",perl=TRUE),"[",1), #regex to split on last _ (sometimes phenotypes have an '_')
                         DAP =       as.numeric(sapply(strsplit(sapply(strsplit(sapply(strsplit(pvalCols,"_\\s*(?=[^_]+$)",perl=TRUE),"[",2),"\\."),"[",2),"-"),"[",1)),
                         Treatment = sapply(strsplit(sapply(strsplit(pvalCols,"_\\s*(?=[^_]+$)",perl=TRUE),"[",1),"_\\s*(?=[^_]+$)",perl=TRUE),"[",2),
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
####Rather than filter by stringent p-value, keep top 50 (with a p < final thresh) from each phenotype/dap/treat###
nKeep <- 50
filteredSNP <- summariseMelt %>% group_by(variable) %>% top_n(n=nKeep,wt=fisherLogP)
filteredSNP <- filteredSNP[filteredSNP$fisherLogP >= -log10(finalThresh),]

filteredSNP$Phenotype = sapply(strsplit(sapply(strsplit(sapply(strsplit(as.character(filteredSNP$variable),"_\\s*(?=[^_]+$)",perl=TRUE),"[",1),"_\\s*(?=[^_]+$)",perl=TRUE),"[",1),"_\\s*(?=[^_]+$)",perl=TRUE),"[",1) #regex to split on last _ (sometimes phenotypes have an '_')
filteredSNP$DAP = as.numeric(sapply(strsplit(sapply(strsplit(as.character(filteredSNP$variable),"_\\s*(?=[^_]+$)",perl=TRUE),"[",1),"_\\s*(?=[^_]+$)",perl=TRUE),"[",2))
filteredSNP$Treatment = sapply(strsplit(sapply(strsplit(sapply(strsplit(as.character(filteredSNP$variable),"_\\s*(?=[^_]+$)",perl=TRUE),"[",1),"_\\s*(?=[^_]+$)",perl=TRUE),"[",1),"_\\s*(?=[^_]+$)",perl=TRUE),"[",2)
filteredSNP$clumpCenter = round((filteredSNP$minBP+filteredSNP$maxBP)/2)
#filteredSNP$Experiment = "IB003toIB004"
filteredSNP$Experiment = "IB005toIB008"

#IB003 <- filteredSNP

finalSNPListAllExps <- rbindlist(list(IB003,filteredSNP))

####Change DAP for single trait measurments to 1 more than the maximum DAP####
single.measure.traits <- c("area.slope",
                           "CR_num",
                           "CR.Angle",
                           "d13C",
                           "CN_ratio",
                           "gN.m2",
                           "gC.m2",
                           "leaf.area",
                           "leaf.weight",
                           "spec.leaf.area",
                           "stom.density")
finalSNPListAllExps$DAP[finalSNPListAllExps$Phenotype %in% single.measure.traits] <- max(finalSNPListAllExps$DAP)+1

#write.table(filterRes,file="../GWASresults/12.clumpColumn.GWASresults.fisherCombinedPval.LTtraits.IB005toIB008.csv",sep=",",col.names=TRUE,row.names=FALSE)
write.table(finalSNPListAllExps,"../results/12.clumpColumn.GWASresults.filtered.ZbrowseReady.ALLtraitsALLexps.csv",sep=",",col.names=TRUE,row.names=FALSE)
