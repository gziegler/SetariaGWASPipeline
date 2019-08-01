rm(list = ls())
if(!require(ionomicsUtils,quietly = TRUE)){
  library(devtools)
  install_github("gziegler/ionomicsUtils")
}
library(ionomicsUtils)
library(data.table)
library(mlmm)
options(scipen = 999)

load("../data/genotype/6.AstleBalding.synbreed.kinship.rda") #loads kinship
load("../data/genotype/6.Eigenstrat.population.structure.50PCs.rda") #loads structData
load("../data/genotype/5.filteredSNPs.noHighCorSNPs.2kbDistThresh.0.5neighborLD.0.975LDfilter.rda")#loads filterHighGeno, neighbors, filterHighInfo, filterHighResults


#phenotype <- read.table("../GWASdatasets/1.Setaria_IR_2016_datsetset_GWAS.BLUPsandBLUEs.csv",sep=",",header=TRUE,stringsAsFactors = FALSE)
phenotype <- read.table("../data/phenotype/9.StomatalDensity.phenotypes.csv",sep=",",header=TRUE,stringsAsFactors = FALSE)
traits <- colnames(phenotype)[2:ncol(phenotype)]

#######Open csv containing info about lines in Genotype file#####
genoInfo <- read.table("../data/genotype/Setaria_597_diversity_samples.csv",sep=",",header=TRUE,stringsAsFactors = FALSE,comment.char = "")
genoInfo$Genotype <- gsub("_setaria_12","",genoInfo$New_name)
length(intersect(genoInfo$Genotype,phenotype$Genotype))
setdiff(phenotype$Genotype,genoInfo$Genotype)
####Just 2 (3 including b100) genotypes not found in genotype file
##writing out table of names to keep for vcf filtering
keepLines <- data.frame(from=sapply(genoInfo$LIB[genoInfo$Genotype %in% c("A10.1",intersect(genoInfo$Genotype,phenotype$Genotype))],function(x){paste(rep(x,4),collapse = "_")}),
                        to=genoInfo$Genotype[genoInfo$Genotype %in% c("A10.1",intersect(genoInfo$Genotype,phenotype$Genotype))],stringsAsFactors = FALSE)
keepLines$to[keepLines$to=="A10.1"] <- "A10"

filterHighGeno <- filterHighGeno[,keepLines$from]
stopifnot(identical(keepLines$from,colnames(filterHighGeno)))
colnames(filterHighGeno) <- keepLines$to

kinship <- kinship[keepLines$from,keepLines$from]
stopifnot(identical(colnames(kinship),keepLines$from))
colnames(kinship) <- keepLines$to
row.names(kinship) <- keepLines$to

structData <- structData[keepLines$from,]
stopifnot(identical(row.names(structData),keepLines$from))
row.names(structData) <- keepLines$to

phenotype <- phenotype[phenotype$Genotype %in% keepLines$to,]

####After filtering only phenotyped lines, refilter by minor allele frequency####
maf <- apply(filterHighGeno,1,function(x) {(length(x[which(x==2)])+0.5*length(x[which(x==1)]))/length(x[!is.na(x)])})
filterHighGeno <- filterHighGeno[which(maf >= 0.05 & maf <= 0.95),]
filterHighInfo <- filterHighInfo[which(maf >= 0.05 & maf <= 0.95),]
#surprisingly this removes very few SNPs, probably because it was already at 0.1 for the whole pop


dir.create("../mlmmTemp", showWarnings = FALSE)
dir.create("../GWASresults", showWarnings = FALSE)
for(i in sample(traits)){
  if(file.exists(paste0("../mlmmTemp/workingon",i))){
    next;
  }else{
    write.table(i,paste0("../mlmmTemp/workingon",i))  
    message("Peforming MLMM on ",i,"\n")
    #Make vector of phenotypes and remove missing
    thisPheno = as.numeric(as.matrix(phenotype[,i]))        
    names(thisPheno) <- phenotype$Genotype
    thisPheno <- thisPheno[which(!(is.na(thisPheno)))]
    commonLines <- intersect(colnames(filterHighGeno),names(thisPheno))
    thisPheno <- thisPheno[names(thisPheno) %in% commonLines]
    this.narrow.Genotype <- filterHighGeno[,colnames(filterHighGeno) %in% commonLines]
    this.narrow.kinship <- kinship[row.names(kinship) %in% commonLines,colnames(kinship) %in% commonLines]
    this.narrow.structure <- structData[row.names(structData) %in% commonLines,]
    
    ####Reorder lines to alphabetical, kinship already is#####
    thisPheno <- thisPheno[order(names(thisPheno))]
    this.narrow.Genotype <- this.narrow.Genotype[,order(colnames(this.narrow.Genotype))]
    row.names(this.narrow.Genotype) <- paste(filterHighInfo$chr,filterHighInfo$pos,sep="_")
    this.narrow.kinship <- this.narrow.kinship[order(row.names(this.narrow.kinship)),order(colnames(this.narrow.kinship))]
    this.narrow.structure <- this.narrow.structure[order(row.names(this.narrow.structure)),]
    
    stopifnot(identical(row.names(this.narrow.kinship),names(thisPheno)))
    stopifnot(identical(colnames(this.narrow.Genotype),names(thisPheno)))
    stopifnot(identical(row.names(this.narrow.structure),names(thisPheno)))
    
    pcVsPheno <- data.frame(PC=numeric(0),pVal=numeric(0),r2=numeric(0))
    pcVsPheno <- rbindlist(list(pcVsPheno,
                                data.frame(1, summary(lm(thisPheno ~ this.narrow.structure[,1]))$coefficients[2,4], summary(lm(thisPheno ~ this.narrow.structure[,1]))$r.squared)))
    #for(pcCol in 2:ncol(this.narrow.structure)){
    for(pcCol in 2:10){
      pcVsPheno <- rbindlist(list(pcVsPheno,
                                  data.frame(pcCol, summary(lm(thisPheno ~ this.narrow.structure[,1:pcCol]))$coefficients[(pcCol+1),4], 
                                             summary(lm(thisPheno ~ this.narrow.structure[,1:pcCol]))$r.squared - 
                                               summary(lm(thisPheno ~ this.narrow.structure[,1:(pcCol-1)]))$r.squared)))
    }
    pcsToUse <- pcVsPheno$PC[which(pcVsPheno$pVal <= 0.001)]
    traitInfo <- data.frame(phenotype=i,genoSize=nrow(this.narrow.Genotype),phenoSize=length(thisPheno),PCs=paste(pcsToUse,collapse=":"))
    print(traitInfo)
    if(length(pcsToUse)==0){
      mygwas <- mlmm(Y=thisPheno,X=t(this.narrow.Genotype),K=this.narrow.kinship,maxsteps = 40,nbchunks = 5)
      traitInfo$PCs <- NA
    }else{
      mygwas <- mlmm_cof(Y=thisPheno,X=t(this.narrow.Genotype),K=this.narrow.kinship,cofs=this.narrow.structure[,pcsToUse],maxsteps = 40,nbchunks = 5)
    }
    
    maxCofModel <- mygwas$pval_step[[length(mygwas$pval_step)]]$out[order(mygwas$pval_step[[length(mygwas$pval_step)]]$out$pval),][1:500,]
    maxCofModel$cofactor[maxCofModel$SNP %in% mygwas$pval_step[[length(mygwas$pval_step)]]$cof] <- 1
    maxCofModel$cofactor[is.na(maxCofModel$cofactor)]<-0  
    cofOrder <- data.frame(cof=mygwas$pval_step[[length(mygwas$pval_step)]]$cof,order=1:length(mygwas$pval_step[[length(mygwas$pval_step)]]$cof))
    outAll <- merge(maxCofModel,cofOrder,by.x="SNP",by.y="cof",all.x=T,all.y=F)
    thisNull <- mygwas$pval_step[[1]]$out[mygwas$pval_step[[1]]$out$SNP %in% outAll$SNP,]
    colnames(thisNull)[2] <- "nullPval"
    outAll <- merge(outAll,thisNull,all.x=T,all.y=F,by="SNP")
    outAll$modelAddedPval <- NA
    for(thisOrd in 1:max(outAll$order,na.rm=T)){
      outAll$modelAddedPval[which(outAll$order==thisOrd)] <- mygwas$pval_step[[thisOrd]]$out$pval[mygwas$pval_step[[thisOrd]]$out$SNP==outAll$SNP[which(outAll$order==thisOrd)]]
    }
    rm(thisNull)
    
    mbonf=subset(mygwas$opt_mbonf$out,pval<= 1e-4)
    mbonf$cofactor[mbonf$SNP %in% mygwas$opt_mbonf$cof]<-1
    mbonf$cofactor[is.na(mbonf$cofactor)]<-0  
    outMbonf <- merge(mbonf,cofOrder,by.x="SNP",by.y="cof",all.x=T,all.y=F)
    outMbonf$order[which(outMbonf$cofactor==0)] <- NA
    
    thisNull <- mygwas$pval_step[[1]]$out[mygwas$pval_step[[1]]$out$SNP %in% outMbonf$SNP,]
    colnames(thisNull)[2] <- "nullPval"
    outMbonf <- merge(outMbonf,thisNull,all.x=T,all.y=F,by="SNP")
    
    if(nrow(outMbonf)>0){
      outMbonf$modelAddedPval <- NA
      if(!all(is.na(outMbonf$order))){
        for(thisOrd in 1:max(outMbonf$order,na.rm=T)){
          outMbonf$modelAddedPval[which(outMbonf$order==thisOrd)] <- mygwas$pval_step[[thisOrd]]$out$pval[mygwas$pval_step[[thisOrd]]$out$SNP==outMbonf$SNP[which(outMbonf$order==thisOrd)]]
        }
      }
    }
    rm(thisNull)
    
    
    
    extBIC=subset(mygwas$opt_extBIC$out,pval<= 1e-4)
    extBIC$cofactor[extBIC$SNP %in% mygwas$opt_extBIC$cof]<-1
    extBIC$cofactor[is.na(extBIC$cofactor)]<-0  
    outExtBIC <- merge(extBIC,cofOrder,by.x="SNP",by.y="cof",all.x=T,all.y=F)
    outExtBIC$order[which(outExtBIC$cofactor==0)] <- NA
    
    thisNull <- mygwas$pval_step[[1]]$out[mygwas$pval_step[[1]]$out$SNP %in% outExtBIC$SNP,]
    colnames(thisNull)[2] <- "nullPval"
    outExtBIC <- merge(outExtBIC,thisNull,all.x=T,all.y=F,by="SNP")
    
    if(nrow(outExtBIC)>0){
      outExtBIC$modelAddedPval <- NA
      if(!all(is.na(outExtBIC$order))){
        for(thisOrd in 1:max(outExtBIC$order,na.rm=T)){
          outExtBIC$modelAddedPval[which(outExtBIC$order==thisOrd)] <- mygwas$pval_step[[thisOrd]]$out$pval[mygwas$pval_step[[thisOrd]]$out$SNP==outExtBIC$SNP[which(outExtBIC$order==thisOrd)]]
        }
      }
    }
    rm(thisNull)
    
    
    nullMod <- subset(mygwas$pval_step[[1]]$out,pval<= 1e-4)
    if(nrow(nullMod)==0){
      nullMod <- cbind(nullMod,data.frame(cofactor=numeric(0),order=numeric(0),nullPval=numeric(0),modelAddedPval=numeric(0),model=character(0)))
    }else{
      nullMod$cofactor <- 0
      nullMod$order <- NA
      nullMod$nullPval <- nullMod$pval
      nullMod$modelAddedPval <- NA
      nullMod$model <- "Null"
    }
    if(nrow(outMbonf)==0){
      outMbonf <- cbind(outMbonf,data.frame(modelAddedPval=numeric(0),model=character(0)))
    }else{
      outMbonf$model <- "Mbonf"
    }  
    if(nrow(outExtBIC)==0){
      outExtBIC <- cbind(outExtBIC,data.frame(modelAddedPval=numeric(0),model=character(0)))
    }else{
      outExtBIC$model <- "ExtBIC"
    } 
    if(nrow(outAll)==0){
      outAll <- cbind(outAll,data.frame(modelAddedPval=numeric(0),model=character(0)))
    }else{
      outAll$model <- "MaxCof"
    }  
    
    everything <- rbind(outAll,outMbonf,outExtBIC,nullMod)
    if(nrow(everything)==0){
      everything <- cbind(everything,data.frame(trait=character(0),nSNPs=numeric(0),nLines=numeric(0)))
    }else{
      everything$trait <- traitInfo$phenotype
      everything$nSNPs <- traitInfo$genoSize
      everything$nLines <- traitInfo$phenoSize
      everything$PCs <- traitInfo$PCs
    }  
    everything <- everything[everything$cofactor==1,]
    write.table(everything,paste0("../GWASresults/10.mlmmResults.",i,".csv"),row.names=FALSE,col.names=TRUE,sep=",")
    rm(mygwas,genoMat,thisKin,thisChunkFiles,traitInfo,thisPheno,outExtBIC,outMbonf,outAll,everything)      
    
  }
}
