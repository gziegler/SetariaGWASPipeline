rm(list = ls())
if(!require(ionomicsUtils,quietly = TRUE)){
  library(devtools)
  install_github("gziegler/ionomicsUtils")
}
# if(!require(lme4qtl,quietly = TRUE)){
#   library(devtools)
#   install_github("variani/lme4qtl")
# }
library(reshape2)
library(ionomicsUtils)
library(data.table)
library(mlmm.gwas)
options(scipen = 999)
keepThresh <- 1 #keep all SNPs less than this pvalue
effectThresh <- 1e-4 #calculate effect size for all SNPs less than this pvalue
PCcofs <- 3 #use top N PCs 
###Code is also present, but not implemented, to use only PCs significantly correlated with trait

load("../data/genotype/6.AstleBalding.synbreed.kinship.rda") #loads kinship
load("../data/genotype/6.Eigenstrat.population.structure.50PCs.rda") #loads structData
load("../data/genotype/5.filteredSNPs.noHighCorSNPs.2kbDistThresh.0.5neighborLD.0.975LDfilter.rda")#loads filterHighGeno, neighbors, filterHighInfo, filterHighResults
rm(neighbors)

#phenotype1 <- read.table("../data/phenotype/2.Setaria_IR_2016_datsetset_GWAS.BLUPsandBLUEs.csv",sep=",",header=TRUE,stringsAsFactors = FALSE)
#phenotype <- read.table("../data/phenotype/0304_setaria_exp.ALL_days.BLUPS.fromCharles.csv",sep=",",header=TRUE,stringsAsFactors = FALSE)
#phenotype <- read.table("../data/phenotype/050708_setaria_exp.ALL_days.BLUPS.fromCharles.csv",sep=",",header=TRUE,stringsAsFactors = FALSE)
phenotype <- read.table("../data/phenotype/all_setaria_medians_for_gwas.fromCharles.csv",sep=",",header=TRUE,stringsAsFactors = FALSE)
####Do just DAP18 and height in wet
phenotype <- phenotype[,c("genotype",grep("height",colnames(phenotype),value=TRUE))]
phenotype <- phenotype[,c("genotype",grep("wet",colnames(phenotype),value=TRUE))]
phenotype <- phenotype[,c("genotype",grep("18",colnames(phenotype),value=TRUE))]
#get to 1 column per treatment/trait (this was for the BLUPs files from Charles)
# meltPhenotype <- reshape2::melt(phenotype,id.vars=c("genotype","trt.day"))
# meltPhenotype <- meltPhenotype[which(!(is.na(meltPhenotype$value))),]
# phenotype <- reshape2::dcast(meltPhenotype,...~trt.day+variable)
colnames(phenotype)[1] <- "Genotype"
#phenotype2 <- read.table("../data/phenotype/9.StomatalDensity.phenotypes.csv",sep=",",header=TRUE,stringsAsFactors = FALSE)
#phenotype <- merge(phenotype1,phenotype2,by="Genotype",all=T)
#phenotype$Awning_dry <- NULL
#phenotype$Awning_wet <- NULL
#phenotype <- read.table("../data/phenotype/rawPhenotypeDatasets/d13C_2018_tbxxxx_v8.csv",sep=",",header=TRUE,stringsAsFactors = FALSE)
#phenotype$X <- NULL
#colnames(phenotype)[2] <- "Genotype"
#phenotype$Genotype <- gsub("TB_","TB_setaria_12_",phenotype$Genotype)
traits <- colnames(phenotype)[2:ncol(phenotype)]

#######Open csv containing info about lines in Genotype file#####
genoInfo <- read.table("../data/genotype/Setaria_597_diversity_samples.csv",sep=",",header=TRUE,stringsAsFactors = FALSE,comment.char = "")
#genoInfo$Genotype <- gsub("_setaria_12","",genoInfo$New_name)
genoInfo$Genotype <- genoInfo$New_name

genoInfo$repLib <- sapply(genoInfo$LIB,function(x){paste(rep(x,4),collapse = "_")})
genoInfo <- genoInfo[which(genoInfo$repLib %in% intersect(genoInfo$repLib,colnames(filterHighGeno))),]

keepLines <- data.frame(from=sapply(genoInfo$LIB[genoInfo$Genotype %in% c("A10.1",intersect(genoInfo$Genotype,phenotype$Genotype))],function(x){paste(rep(x,4),collapse = "_")}),
                        to=genoInfo$Genotype[genoInfo$Genotype %in% c("A10.1",intersect(genoInfo$Genotype,phenotype$Genotype))],stringsAsFactors = FALSE)
#keepLines$to[keepLines$to=="A10.1"] <- "A10"

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

####Add an X to the beginning of each SNP id, one of the functions below doesn't like it starting with a digit####
row.names(filterHighGeno) <- paste0("X",filterHighInfo$chr,"_",filterHighInfo$pos)
rm(filterHighInfo)

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
    if(length(commonLines)<10){
      warning("Phenotype ",i," only has ",length(commonLines)," genotyped observations. Not running GWAS")
      next;
    }
    thisPheno <- thisPheno[names(thisPheno) %in% commonLines]
    this.narrow.Genotype <- filterHighGeno[,colnames(filterHighGeno) %in% commonLines]
    this.narrow.kinship <- kinship[row.names(kinship) %in% commonLines,colnames(kinship) %in% commonLines]
    this.narrow.structure <- structData[row.names(structData) %in% commonLines,]
    
    ####Reorder lines to alphabetical, kinship already is#####
    thisPheno <- thisPheno[order(names(thisPheno))]
    this.narrow.Genotype <- this.narrow.Genotype[,order(colnames(this.narrow.Genotype))]
    this.narrow.kinship <- this.narrow.kinship[order(row.names(this.narrow.kinship)),order(colnames(this.narrow.kinship))]
    this.narrow.structure <- this.narrow.structure[order(row.names(this.narrow.structure)),]
    
    stopifnot(identical(row.names(this.narrow.kinship),names(thisPheno)))
    stopifnot(identical(colnames(this.narrow.Genotype),names(thisPheno)))
    stopifnot(identical(row.names(this.narrow.structure),names(thisPheno)))
    
    pcVsPheno <- data.frame(PC=numeric(0),pVal=numeric(0),r2=numeric(0))
    pcVsPheno <- rbindlist(list(pcVsPheno,
                                data.frame(PC=1, summary(lm(thisPheno ~ this.narrow.structure[,1]))$coefficients[2,4], summary(lm(thisPheno ~ this.narrow.structure[,1]))$r.squared)),use.names=FALSE)
    #for(pcCol in 2:ncol(this.narrow.structure)){
    for(pcCol in 2:10){
      pcVsPheno <- rbindlist(list(pcVsPheno,
                                  data.frame(pcCol, summary(lm(thisPheno ~ this.narrow.structure[,1:pcCol]))$coefficients[(pcCol+1),4], 
                                             summary(lm(thisPheno ~ this.narrow.structure[,1:pcCol]))$r.squared - 
                                               summary(lm(thisPheno ~ this.narrow.structure[,1:(pcCol-1)]))$r.squared)),use.names = F)
    }
    pcsToUse <- pcVsPheno$PC[which(pcVsPheno$pVal <= 0.001)]
    #traitInfo <- data.frame(phenotype=i,genoSize=nrow(this.narrow.Genotype),phenoSize=length(thisPheno),PCs=paste(pcsToUse,collapse=":"))
    traitInfo <- data.frame(phenotype=i,genoSize=nrow(this.narrow.Genotype),phenoSize=length(thisPheno),PCs=paste(1:PCcofs,collapse=":"))
    print(traitInfo)
    # if(length(pcsToUse)==0){
    #   mygwas <- mlmm(Y=thisPheno,X=t(this.narrow.Genotype),K=this.narrow.kinship,maxsteps = 40,nbchunks = 5)
    #   traitInfo$PCs <- NA
    # }else{
    #   mygwas <- mlmm_cof(Y=thisPheno,X=t(this.narrow.Genotype),K=this.narrow.kinship,cofs=this.narrow.structure[,pcsToUse],maxsteps = 40,nbchunks = 5)
    # }
    # 
    colnames(this.narrow.structure) <- paste0("PC",c(1:ncol(this.narrow.structure)))
    this.narrow.structure <- this.narrow.structure[,1:PCcofs]
    
    #mygwas <- mlmm(Y=thisPheno,X=t(this.narrow.Genotype),K=this.narrow.kinship,maxsteps = 3,nbchunks = 3)
    #mygwas <- mlmm_cof(Y=thisPheno,X=t(this.narrow.Genotype),K=this.narrow.kinship,cofs=this.narrow.structure[,1:3],maxsteps = 3,nbchunks = 3)
    
    mygwasNew <- mlmm_allmodels(Y=thisPheno,XX=list(t(this.narrow.Genotype)),KK=list(this.narrow.kinship),cofs=this.narrow.structure,maxsteps = 3,nbchunks = 3)
    #mygwasnocof <- mlmm_allmodels(Y=thisPheno,XX=list(t(this.narrow.Genotype)),KK=list(this.narrow.kinship),maxsteps = 3,nbchunks = 3)
    
    #####Pipeline for getting effect sizes using EBIC#####
    #pulls out SNPs from genotype file that were selected to be included in the model)
    #sel_XX <- frommlmm_toebic(list(t(this.narrow.Genotype)),mygwasNew)
    #
    #calculate ebic for all models
    #res.eBIC <- eBIC_allmodels(thisPheno, sel_XX, list(this.narrow.kinship), nb.tests = nrow(this.narrow.Genotype),cofs = this.narrow.structure[,1:3])
    # 
    #res.threshold <- threshold_allmodels(threshold = 0.05,mygwasNew)
    # 
    #sel_XXclass <- fromeBICtoEstimation(list(t(this.narrow.Genotype)), res.eBIC,res.threshold=res.threshold)
    #sel_XXclass$X2_27572311 <- NULL
    # 
    #eff.estimations <- Estimation_allmodels(thisPheno,sel_XXclass,list(this.narrow.kinship),this.narrow.structure[,1:3])
    #genotypes.boxplot(t(this.narrow.Genotype),thisPheno,effects=eff.estimations)
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #Code to get effect size for each individual SNP above some threshold################
    #####################################################################################
    #####################################################################################
    ####Gets pvals from the single SNP model (regular MLM)#### 
    keepSNPs <- mygwasNew[2][[1]][which(mygwasNew[2][[1]]<=keepThresh)]
    keepSNPs <- data.frame(SNP=names(keepSNPs),pval=keepSNPs,cofactor=0,model="Single",
                           trait=traitInfo$phenotype,nSNPs=traitInfo$genoSize,nLines=traitInfo$phenoSize,PCs=traitInfo$PCs,stringsAsFactors = FALSE)
    keepSNPs$effectSize <- NA #column to fill with the effect sizes calculated below
    keepSNPs$maf <- NA
    
    #####Use these lines to calculate pseudoH####
    ######Calculate pseudo-H of model with just kinship and cofactors#######
    #n <- length(thisPheno)
    # X0<-as.matrix(rep(1,n))
    # colnames(X0)<-"mu"
    # #cof <- cbind(X0,this.narrow.structure[,1:3],data.frame("7_3771623"=t.this.narrow.Genotype[,"7_3771623"],"7_3337761"=t.this.narrow.Genotype[,"7_3337761"]))
    # cof <- cbind(X0,this.narrow.structure[,1:3])
    # ind<-as.factor(names(thisPheno)) #default for additive and dominance model
    # 
    # data = data.frame(id=names(thisPheno), Y=as.vector(thisPheno), ind, cof)
    # fixed = formula(sub("mu", "1", paste("Y~",paste(colnames(cof), collapse = "+"))))
    # 
    # #kinship normalization - orig mlmm did this too
    #cst<-(n-1)/sum((diag(n)-matrix(1,n,n)/n)*this.narrow.kinship )
    #KK.norm <- cst*this.narrow.kinship
    
    #There is a step in mlmm to make matrix positive semi-definite, but I think the way this one is calculated it should always be
    # res = sommer::mmer(fixed=fixed, random=~vs(id, Gu=KK.norm), data=data) #can add method = "emma" to get h2 almost exactly as they are in orig MLMM
    # h2 <- 1-as.vector(res$sigma[['units']])/Reduce("+", res$sigma) #It's estimating H2 as 1-kinV/totalV
    # keepSNPs$pseudoh <- as.numeric(h2)
    effectSNPs <- which(keepSNPs$pval <= effectThresh)
    for(snp in effectSNPs){
      message("Calculating effect size for ",keepSNPs$SNP[snp])
      thisSNP <- as.data.frame(this.narrow.Genotype[keepSNPs$SNP[snp],])
      colnames(thisSNP) <- keepSNPs$SNP[snp]
      thisSNP[thisSNP[,1]==0,1] <- "00"
      thisSNP[thisSNP[,1]==1,1] <- "01|10"
      thisSNP[thisSNP[,1]==2,1] <- "11"
      thisSNP[,1] <- factor(thisSNP[,1])
      eff.estimations <- Estimation_allmodels(thisPheno,thisSNP,list(this.narrow.kinship),this.narrow.structure)
      #microbenchmark(effEst = {eff.estimations <- Estimation_allmodels(thisPheno,thisSNP,list(this.narrow.kinship),this.narrow.structure)},
      #               time=10,unit="relative")
      ###Using lme4qtl gives the same effect sizes
      #design <- data.frame(Y=thisPheno,X0,this.narrow.structure[,1:3],thisSNP,ind=as.factor(names(thisPheno)))
      #lmerRes <- lme4qtl::relmatLmer(Y ~ 0 + mu + PC1 + PC2 + PC3 + X7_19588148 + (1|ind),design,relmat=list(ind=this.narrow.kinship))
      
      SNPeffect <- eff.estimations[paste0(keepSNPs$SNP[snp],"_00"),"BLUE"]*-1
      if(SNPeffect==0){
        warning("Minor allele wasn't the mean effect")
        SNPeffect <- eff.estimations[paste0(keepSNPs$SNP[snp],"_11"),"BLUE"]
      }
      keepSNPs$effectSize[snp] <- SNPeffect
      keepSNPs$maf[snp] <- eff.estimations[paste0(keepSNPs$SNP[snp],"_11"),"Frequency"]
      # #####below is code that is what is implemented in Estimation_allmodels#####
      # #cofs is this.narrow.structure[,1:3]
      # n <- length(thisPheno)
      # X0<-as.matrix(rep(1,n))
      # colnames(X0)<-"mu"
      # ind<-as.factor(names(thisPheno)) #default for additive and dominance model
      # 
      # #cof <- cbind(X0,this.narrow.structure[,1:3],data.frame("7_3771623"=t.this.narrow.Genotype[,"7_3771623"],"7_3337761"=t.this.narrow.Genotype[,"7_3337761"]))
      # #design matrix (data frame of all fixed effects)
      # design <- data.frame(X0,this.narrow.structure[,1:3],thisSNP)
      # #design[,which(!(colnames(design) %in% c("mu",colnames(this.narrow.structure))))]
      # randForm <- c("(1|ind)")
      # xm <- as.formula(paste(c("Y~0",colnames(design),randForm),sep="+",collapse = "+"))
      # 
      # res.lmekin<-coxme::lmekin(xm, data=data.frame(Y=thisPheno, design, ind),
      #                           varlist = coxme::coxmeMlist(KK.norm) ,method="ML",
      #                           model = TRUE, x = TRUE, y = TRUE,
      #                           na.action = na.omit)
      # blue = res.lmekin$coefficients$fixed
      ####This point in Estimation_allmodels_withlmekin checks for non-estimated levels in the fixed part of the model####
      ##If that becomes a problem then that section should be added in, starts at line 174 of that code
      
    }#end calculate effect sizes for SNPs
    keepSNPs$SNP <- sub("^X","",keepSNPs$SNP)
    #write.table(keepSNPs,paste0("../GWASresults/10.mlmmResults-",i,".csv"),row.names=FALSE,col.names=TRUE,sep=",")
    ##Write in highly compressed rda, write_parquet is crazy fast, but needs to be installed on the server
    save(keepSNPs,file=paste0("../GWASresults/10.mlmmResults-",i,".rda"),compress="xz",compression_level = 5)
    #arrow::write_parquet(keepSNPs,sink=paste0("../GWASresults/10.mlmmResults-",i,".parquet"),compression="gzip",compression_level = 5)
    #parqSNPs <- arrow::read_parquet(file=paste0("../GWASresults/10.mlmmResults-",i,".parquet"),as_data_frame = TRUE)
    rm(mygwasNew,keepSNPs,this.narrow.Genotype,thisPheno,this.narrow.kinship,this.narrow.structure)      
    
    
    # maxCofModel <- mygwas$pval_step[[length(mygwas$pval_step)]]$out[order(mygwas$pval_step[[length(mygwas$pval_step)]]$out$pval),][1:5000,]
    # maxCofModel$cofactor[maxCofModel$SNP %in% mygwas$pval_step[[length(mygwas$pval_step)]]$cof] <- 1
    # maxCofModel$cofactor[is.na(maxCofModel$cofactor)]<-0  
    # cofOrder <- data.frame(cof=mygwas$pval_step[[length(mygwas$pval_step)]]$cof,order=1:length(mygwas$pval_step[[length(mygwas$pval_step)]]$cof))
    # outAll <- merge(maxCofModel,cofOrder,by.x="SNP",by.y="cof",all.x=T,all.y=F)
    # thisNull <- mygwas$pval_step[[1]]$out[mygwas$pval_step[[1]]$out$SNP %in% outAll$SNP,]
    # colnames(thisNull)[2] <- "nullPval"
    # outAll <- merge(outAll,thisNull,all.x=T,all.y=F,by="SNP")
    # outAll$modelAddedPval <- NA
    # for(thisOrd in 1:max(outAll$order,na.rm=T)){
    #   outAll$modelAddedPval[which(outAll$order==thisOrd)] <- mygwas$pval_step[[thisOrd]]$out$pval[mygwas$pval_step[[thisOrd]]$out$SNP==outAll$SNP[which(outAll$order==thisOrd)]]
    # }
    # rm(thisNull)
    # 
    # mbonf=subset(mygwas$opt_mbonf$out,pval<= 1e-3)
    # mbonf$cofactor[mbonf$SNP %in% mygwas$opt_mbonf$cof]<-1
    # mbonf$cofactor[is.na(mbonf$cofactor)]<-0
    # outMbonf <- merge(mbonf,cofOrder,by.x="SNP",by.y="cof",all.x=T,all.y=F)
    # outMbonf$order[which(outMbonf$cofactor==0)] <- NA
    # 
    # thisNull <- mygwas$pval_step[[1]]$out[mygwas$pval_step[[1]]$out$SNP %in% outMbonf$SNP,]
    # colnames(thisNull)[2] <- "nullPval"
    # outMbonf <- merge(outMbonf,thisNull,all.x=T,all.y=F,by="SNP")
    # 
    # if(nrow(outMbonf)>0){
    #   outMbonf$modelAddedPval <- NA
    #   if(!all(is.na(outMbonf$order))){
    #     for(thisOrd in 1:max(outMbonf$order,na.rm=T)){
    #       outMbonf$modelAddedPval[which(outMbonf$order==thisOrd)] <- mygwas$pval_step[[thisOrd]]$out$pval[mygwas$pval_step[[thisOrd]]$out$SNP==outMbonf$SNP[which(outMbonf$order==thisOrd)]]
    #     }
    #   }
    # }
    # rm(thisNull)
    # 
    # 
    # 
    # extBIC=subset(mygwas$opt_extBIC$out,pval<= 1e-3)
    # extBIC$cofactor[extBIC$SNP %in% mygwas$opt_extBIC$cof]<-1
    # extBIC$cofactor[is.na(extBIC$cofactor)]<-0  
    # outExtBIC <- merge(extBIC,cofOrder,by.x="SNP",by.y="cof",all.x=T,all.y=F)
    # outExtBIC$order[which(outExtBIC$cofactor==0)] <- NA
    # 
    # thisNull <- mygwas$pval_step[[1]]$out[mygwas$pval_step[[1]]$out$SNP %in% outExtBIC$SNP,]
    # colnames(thisNull)[2] <- "nullPval"
    # outExtBIC <- merge(outExtBIC,thisNull,all.x=T,all.y=F,by="SNP")
    # 
    # if(nrow(outExtBIC)>0){
    #   outExtBIC$modelAddedPval <- NA
    #   if(!all(is.na(outExtBIC$order))){
    #     for(thisOrd in 1:max(outExtBIC$order,na.rm=T)){
    #       outExtBIC$modelAddedPval[which(outExtBIC$order==thisOrd)] <- mygwas$pval_step[[thisOrd]]$out$pval[mygwas$pval_step[[thisOrd]]$out$SNP==outExtBIC$SNP[which(outExtBIC$order==thisOrd)]]
    #     }
    #   }
    # }
    # rm(thisNull)
    # 
    # 
    # nullMod <- subset(mygwas$pval_step[[1]]$out,pval<= 1)
    # if(nrow(nullMod)==0){
    #   nullMod <- cbind(nullMod,data.frame(cofactor=numeric(0),order=numeric(0),nullPval=numeric(0),modelAddedPval=numeric(0),model=character(0)))
    # }else{
    #   nullMod$cofactor <- 0
    #   nullMod$order <- NA
    #   nullMod$nullPval <- nullMod$pval
    #   nullMod$modelAddedPval <- NA
    #   nullMod$model <- "Null"
    # }
    # if(nrow(outMbonf)==0){
    #   outMbonf <- cbind(outMbonf,data.frame(modelAddedPval=numeric(0),model=character(0)))
    # }else{
    #   outMbonf$model <- "Mbonf"
    # }  
    # if(nrow(outExtBIC)==0){
    #   outExtBIC <- cbind(outExtBIC,data.frame(modelAddedPval=numeric(0),model=character(0)))
    # }else{
    #   outExtBIC$model <- "ExtBIC"
    # } 
    # if(nrow(outAll)==0){
    #   outAll <- cbind(outAll,data.frame(modelAddedPval=numeric(0),model=character(0)))
    # }else{
    #   outAll$model <- "MaxCof"
    # }  
    # 
    # everything <- rbind(outAll,outMbonf,outExtBIC,nullMod)
    # if(nrow(everything)==0){
    #   everything <- cbind(everything,data.frame(trait=character(0),nSNPs=numeric(0),nLines=numeric(0)))
    # }else{
    #   everything$trait <- traitInfo$phenotype
    #   everything$nSNPs <- traitInfo$genoSize
    #   everything$nLines <- traitInfo$phenoSize
    #   everything$PCs <- traitInfo$PCs
    # }  
    #everything <- everything[everything$cofactor==1,]
    # write.table(everything,paste0("../GWASresults/10.mlmmResults.",i,".csv"),row.names=FALSE,col.names=TRUE,sep=",")
    # rm(mygwas,genoMat,thisKin,thisChunkFiles,traitInfo,thisPheno,outExtBIC,outMbonf,outAll,everything)      
    # 
  }
}
