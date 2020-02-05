rm(list=ls())
library(data.table)
library(dplyr)

#first pass fisherP threshold - filter SNPs by this fisherP, next pass will clump remaining SNPs and select best p-val from clumps
firstThresh <- 1e-2

#library(metap) #for fisher's combined
resultsFiles <- list.files("../GWASresults/","\\.rda",full.names=TRUE)
#resultsFiles <- grep("IB005|IB007|IB008",resultsFiles,value=TRUE)
#resultsFiles <- grep("IB003|IB004",resultsFiles,value=TRUE)
resultsFiles <- grep("rank",resultsFiles,value=TRUE)
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

#Deconvulte Trait ID into a table
traitTable <- data.frame(Phenotype = sapply(strsplit(colnames(mergeRes)[2:ncol(mergeRes)],"_"),"[",1),
                         DAP =       as.numeric(sapply(strsplit(sapply(strsplit(sapply(strsplit(colnames(mergeRes)[2:ncol(mergeRes)],"_"),"[",2),"\\."),"[",3),"-"),"[",1)),
                         Treatment = sapply(strsplit(sapply(strsplit(colnames(mergeRes)[2:ncol(mergeRes)],"_"),"[",2),"\\."),"[",2),
                         Experiment = sapply(strsplit(sapply(strsplit(colnames(mergeRes)[2:ncol(mergeRes)],"_"),"[",2),"\\."),"[",1),
                         colID = colnames(mergeRes)[2:ncol(mergeRes)],
                         stringsAsFactors = FALSE)
#calculate fishers combined on each Phenotype/DAP/treatment combination
#create function for fisher combining pvals, function from Charles Pignon, gives same answer as metap::sumlog
fisher_combine_pvals=function(pval_list){pchisq((sum(log(pval_list))*-2), df=length(pval_list)*2, lower.tail=F)}
mergeRes <- as.data.table(mergeRes)
for(pheno in unique(traitTable$Phenotype)){
  for(dap in unique(traitTable$DAP[which(traitTable$Phenotype==pheno)])){
    for(treat in unique(traitTable$Treatment[which(traitTable$Phenotype==pheno & traitTable$DAP==dap)])){
      message(pheno," ",dap," ",treat,"\n")
      thisTraitCols <- grep("pval",traitTable$colID[which(traitTable$Phenotype==pheno & traitTable$DAP==dap & traitTable$Treatment==treat)],value=TRUE)
      #thisTraitCols <- grep("mean",thisTraitCols,value=TRUE,invert = TRUE) #exclude columns from rank.mean experiment
      mergeRes[,paste(pheno,treat,dap,"fisherP",sep="_") := fisher_combine_pvals(.SD),by=1:nrow(mergeRes),.SDcols = thisTraitCols]
    }
  }
}

fisherPcols <- grep("fisherP",colnames(mergeRes),value=TRUE)
#fisherPcols <- grep("fisherP|mean.*pval",colnames(mergeRes),value=TRUE) #don't filter out mean pval column (for the rank dataset)
filterRes <- mergeRes[mergeRes[, Reduce(`|`, lapply(.SD, `<=`, firstThresh)),.SDcols = fisherPcols]] #keep phenotype columns where fisherP <= firstThresh
#save(mergeRes,file=paste0("../GWASresults/11.CombinedMLMMResultsWithFisherPval.allSNPs.rda"),compress="xz",compression_level = 5)
#save(filterRes,file=paste0("../GWASresults/11.CombinedMLMMResultsWithFisherPval.filtered.rda"),compress="xz",compression_level = 5)

#save(mergeRes,file=paste0("../GWASresults/11.JustIB005to8.CombinedMLMMResultsWithFisherPval.allSNPs.rda"),compress="xz",compression_level = 5)
#save(filterRes,file=paste0("../GWASresults/11.JustIB005to8.CombinedMLMMResultsWithFisherPval.filtered.rda"),compress="xz",compression_level = 5)
# save(mergeRes,file=paste0("../GWASresults/11.JustIB003to4.CombinedMLMMResultsWithFisherPval.allSNPs.rda"),compress="xz",compression_level = 5)
# save(filterRes,file=paste0("../GWASresults/11.JustIB003to4.CombinedMLMMResultsWithFisherPval.filtered.rda"),compress="xz",compression_level = 5)
save(mergeRes,file=paste0("../GWASresults/11.RANK5to8.CombinedMLMMResultsWithFisherPval.allSNPs.rda"),compress="xz",compression_level = 5)
save(filterRes,file=paste0("../GWASresults/11.RANK5to8.CombinedMLMMResultsWithFisherPval.filtered.rda"),compress="xz",compression_level = 5)

####Move on to clump SNPs####





#recalculate combined by taking best pval from each clump foreach experiment...

#do effect size calculation on high pval SNPs to filter further...