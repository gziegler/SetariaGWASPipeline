rm(list=ls())
library(data.table)
library(dplyr)
library(feather)
###To save time this doesn't calcualte fisherP for every SNP, it first filters any SNPs where the minimum pvalue (across all phenotypes) is > 0.1
#####A pvalue of 0.06 across 3 experiments gives a fisher p of 0.0097
######This actually only removes about 5000 SNPs for IB005,IB007,IB008
##So ended up not doing it
#calcFisher <- 0.06 #threshold to filter SNPs before calculating fisher's
#first pass fisherP threshold - filter SNPs by this fisherP, next pass will clump remaining SNPs and select best p-val from clumps
firstThresh <- 1e-2

#library(metap) #for fisher's combined
resultsFiles <- list.files("../GWASresults/","\\.rda",full.names=TRUE)
resultsFiles <- grep("10.mlmmResults",resultsFiles,value=TRUE)
resultsFiles <- grep("IB005|IB007|IB008",resultsFiles,value=TRUE)
#resultsFiles <- grep("IB003|IB004",resultsFiles,value=TRUE)
#resultsFiles <- grep("rank",resultsFiles,value=TRUE)
mergeRes <- data.frame()
for(i in resultsFiles){
  message("Reading results file ",i,"\n")
  load(i)#these contain object keepSNPs
  keepSNPs[,which(!(colnames(keepSNPs) %in% c("SNP","pval","effectSize","trait")))] <- NULL
  setnames(keepSNPs,old=c('pval','effectSize'),new=c(paste(keepSNPs$trait[1],"pval",sep="-"),paste(keepSNPs$trait[1],"effectSize",sep="-")))
  keepSNPs$trait <- NULL
  if(nrow(mergeRes)==0){
    mergeRes <- keepSNPs
  }else{
    if(identical(mergeRes$SNP,keepSNPs$SNP)){
      keepSNPs$SNP <- NULL
      mergeRes <- cbind(mergeRes,keepSNPs)
    }else{
      stop("SNP order not identical between files.")
    }
  }
  rm(keepSNPs)
}

#can write it out here to save time later, same file + fisherP is written out below
#save(mergeRes,file="../GWASresults/11.IB005toIB008mergedResults.rda",compress="xz",compression_level = 5)
#arrow::write_parquet(mergeRes,sink="../GWASresults/11.IB005toIB008mergedResults.parquet",compression="gzip",compression_level = 5)
#feather::write_feather(mergeRes,path="../GWASresults/11.IB005toIB008mergedResults.feather")
#arrow::write_parquet(mergeRes,sink="../GWASresults/11.IB003toIB004mergedResults.parquet",compression="gzip",compression_level = 5)
#feather::write_feather(mergeRes,path="../GWASresults/11.IB003toIB004mergedResults.feather")

#Deconvulte Trait ID into a table
traitTable <- data.frame(Phenotype = sapply(strsplit(sapply(strsplit(colnames(mergeRes)[2:ncol(mergeRes)],"_\\s*(?=[^_]+$)",perl=TRUE),"[",1),"_\\s*(?=[^_]+$)",perl=TRUE),"[",1), #regex to split on last _ (sometimes phenotypes have an '_')
                         DAP =       as.numeric(sapply(strsplit(sapply(strsplit(sapply(strsplit(colnames(mergeRes)[2:ncol(mergeRes)],"_\\s*(?=[^_]+$)",perl=TRUE),"[",2),"\\."),"[",2),"-"),"[",1)),
                         Treatment = sapply(strsplit(sapply(strsplit(colnames(mergeRes)[2:ncol(mergeRes)],"_\\s*(?=[^_]+$)",perl=TRUE),"[",1),"_\\s*(?=[^_]+$)",perl=TRUE),"[",2),
                         Experiment = sapply(strsplit(sapply(strsplit(colnames(mergeRes)[2:ncol(mergeRes)],"_\\s*(?=[^_]+$)",perl=TRUE),"[",2),"\\."),"[",1),
                         colID = colnames(mergeRes)[2:ncol(mergeRes)],
                         stringsAsFactors = FALSE)
#calculate fishers combined on each Phenotype/DAP/treatment combination
#create function for fisher combining pvals, function from Charles Pignon, gives same answer as metap::sumlog
#fisher_combine_pvals=function(pval_list){pchisq((sum(log(pval_list))*-2), df=length(pval_list)*2, lower.tail=F)} #calculating this in data.table is 
###much, much faster
# message(nrow(mergeRes)," SNPs before first filter.")
# setDT(mergeRes)
# pvalCols <- grep("pval",colnames(mergeRes),value=TRUE)
# mergeRes[,row.min := do.call(pmin, c(.SD,na.rm = TRUE)),.SDcols=pvalCols]
# firstFilter <- mergeRes[row.min < calcFisher] #keep rows where minimum is less than the threshold
# firstFilter$row.min <- NULL
# message(nrow(firstFilter)," SNPs after first filter.")
# mergeRes <- firstFilter
# rm(firstFilter)
#For all of the lemnatec traits, this apparently ended up removing almost no SNPs
setDT(mergeRes)
for(pheno in unique(traitTable$Phenotype)){
  for(dap in unique(traitTable$DAP[which(traitTable$Phenotype==pheno)])){
    for(treat in unique(traitTable$Treatment[which(traitTable$Phenotype==pheno & traitTable$DAP==dap)])){
      message(pheno," ",dap," ",treat,"\n")
      thisTraitCols <- grep("pval",traitTable$colID[which(traitTable$Phenotype==pheno & traitTable$DAP==dap & traitTable$Treatment==treat)],value=TRUE)
      chiSqDF = 2*length(thisTraitCols)
      #thisTraitCols <- grep("mean",thisTraitCols,value=TRUE,invert = TRUE) #exclude columns from rank.mean experiment
      mergeRes[,paste(pheno,treat,dap,"fisherP",sep="_") := pchisq((-2*rowSums(log(.SD),na.rm=T)),df=chiSqDF,lower.tail=F), .SDcols = thisTraitCols]
      #mergeRes[,paste(pheno,treat,dap,"fisherP",sep="_") := fisher_combine_pvals(.SD),by=seq_len(nrow(mergeRes)),.SDcols = thisTraitCols]
    }
  }
}

fisherPcols <- grep("fisherP",colnames(mergeRes),value=TRUE)
#fisherPcols <- grep("fisherP|mean.*pval",colnames(mergeRes),value=TRUE) #don't filter out mean pval column (for the rank dataset)
filterRes <- mergeRes[mergeRes[, Reduce(`|`, lapply(.SD, `<=`, firstThresh)),.SDcols = fisherPcols]] #keep phenotype columns where fisherP <= firstThresh
#save(mergeRes,file=paste0("../GWASresults/11.CombinedMLMMResultsWithFisherPval.allSNPs.rda"),compress="xz",compression_level = 5)
#save(filterRes,file=paste0("../GWASresults/11.allLTphenos.IB005to8.CombinedMLMMResultsWithFisherPval.filtered.rda"),compress="xz",compression_level = 5)
#feather::write_feather(mergeRes,path="../GWASresults/11.allLTphenos.IB003to4.CombinedMLMMResultsWithFisherPval.allSNPs.feather")
feather::write_feather(mergeRes,path="../GWASresults/11.allLTphenos.IB005to8.CombinedMLMMResultsWithFisherPval.allSNPs.feather")

#save(mergeRes,file=paste0("../GWASresults/11.JustIB005to8.CombinedMLMMResultsWithFisherPval.allSNPs.rda"),compress="xz",compression_level = 5)
#save(filterRes,file=paste0("../GWASresults/11.JustIB005to8.CombinedMLMMResultsWithFisherPval.filtered.rda"),compress="xz",compression_level = 5)
# save(mergeRes,file=paste0("../GWASresults/11.JustIB003to4.CombinedMLMMResultsWithFisherPval.allSNPs.rda"),compress="xz",compression_level = 5)
# save(filterRes,file=paste0("../GWASresults/11.JustIB003to4.CombinedMLMMResultsWithFisherPval.filtered.rda"),compress="xz",compression_level = 5)
#save(mergeRes,file=paste0("../GWASresults/11.RANK5to8.CombinedMLMMResultsWithFisherPval.allSNPs.rda"),compress="xz",compression_level = 5)
#save(filterRes,file=paste0("../GWASresults/11.RANK5to8.CombinedMLMMResultsWithFisherPval.filtered.rda"),compress="xz",compression_level = 5)

####Move on to clump SNPs####





#recalculate combined by taking best pval from each clump foreach experiment...

#do effect size calculation on high pval SNPs to filter further...