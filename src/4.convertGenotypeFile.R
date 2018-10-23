rm(list=ls())
library(data.table)
if(!require(ionomicsUtils,quietly = TRUE)){
  library(devtools)
  install_github("gziegler/ionomicsUtils")
}
library(ionomicsUtils)

###There are no missing because this was the imputed set from Sujan
###MAF filtering already happened in vcftools, so there are also no SNPs removed for MAF
hetThreshold <- 0.25 #Because of the density of SNPs, I REMOVE SNPs with a het value greater than this threshold
maf <- 0.1
maxMissingBySNP <- 0.2

geno <- fread(paste0("../data/genotype/3.from12.setaria.maf0.1.maxMissing0.1.012"),header=FALSE,sep="\t")
indv <- read.table(paste0("../data/genotype/3.from12.setaria.maf0.1.maxMissing0.1.012.indv"),sep="\t",header=FALSE,stringsAsFactors = FALSE)
snpInfo <- fread(paste0("../data/genotype/3.from12.setaria.maf0.1.maxMissing0.1.012.pos"),sep="\t",header=FALSE,stringsAsFactors = FALSE)
setnames(snpInfo,c("V1","V2"),c("chr","pos"))
#First column of geno is  rownumbers 0:length(lines)
geno$V1 <- indv$V1
setnames(geno,colnames(geno),c("Genotype",paste0(snpInfo$chr,"_",snpInfo$pos)))

geno <- transpose(geno)
setnames(geno, old = names(geno), new = as.character(geno[1,]))
geno <- geno[-1,]

###Change -1 values to missing
for(j in seq_len(ncol(geno))){
  set(geno,which(geno[[j]]==-1),j,NA)
}

##geno now has columns as individuals, rows are SNPs, there are no metadata columns
alleleTable <- cbind(snpInfo,ionomicsUtils::calcHet(geno))

# alleleTable[,discard := (FracHet>=hetThreshold)]
# if(length(which(alleleTable$discard==TRUE))>0){
#   geno[which(alleleTable$discard==TRUE)][geno[which(alleleTable$discard==TRUE)]==1] <- NA
# }
# alleleTable[, discard := NULL]

###############Filter for minor allele frequency, missing SNPs, etc.###########
alleleTable <- cbind(alleleTable,ionomicsUtils::calcMAF(geno))
#decide who to delete
alleleTable[,discard := (FracHet>hetThreshold | FracMissing>maxMissingBySNP | MAF<maf | MAF>(1-maf))]

#get rid of the rows in snp that are TRUE in alleleTable$discard
keepIdxs <- which(alleleTable$discard==FALSE)
cols <- names(geno)
firstCol <- cols[1]
snp.subset <- data.table(col1 = geno[[firstCol]][keepIdxs]) 
setnames(snp.subset,"col1",firstCol)
for(col in cols[2:length(cols)]){
  snp.subset[, (col) := geno[[col]][keepIdxs],]
  geno[, (col):= NULL,] #delete
}
rm(geno)
alleleTable <- alleleTable[discard==FALSE]
alleleTable[, discard := NULL]

####Impute missing values to major allele###
######Note that in this snp file not all 0s are the major allele####
#######THERE ARE NO MISSING VALUES DUE TO PRIOR IMPUTATION#####
####So, check if MAF >0.5, impute to 2
###Change -1 values to missing
# for(j in seq_len(ncol(geno))){
#   set(geno,which(is.na(geno[[j]]) & alleleTable$MAF<=0.5),j,0)
#   set(geno,which(is.na(geno[[j]]) & alleleTable$MAF>0.5),j,2)
# }
# 


save(snp.subset,alleleTable,file="../data/genotype/4.FilteredGenotypeFile.hetFilter0.25.maf0.1.rda")

#we can now convert to a matrix and use the apply functions
genoMatrix <- as.matrix(snp.subset)
rm(snp.subset)
class(genoMatrix) <- "numeric"

alleleTable[, chr := sapply(strsplit(alleleTable$chr, "_"), "[", 2)]
alleleTable[, chr := as.numeric(chr)]

keepIdxs <- which(alleleTable$chr %in% c(1:9)) #matches first 1:4420660 indices identical(keepIdxs[1:4236621],c(1:4236621))
genoMatrix <- genoMatrix[keepIdxs,]
alleleTable <- alleleTable[keepIdxs,]
rm(keepIdxs)

#Resort snpInfo and genoMatrix so chromsome/bp are in order
keepIdxs <- order(alleleTable$chr,alleleTable$pos)
genoMatrix <- genoMatrix[keepIdxs,]
alleleTable <- alleleTable[keepIdxs,]
rm(keepIdxs)

save(genoMatrix,alleleTable,file="../data/genotype/4.FilteredGenotypeFile.MatrixFormat.noscaffold.hetFilter0.25.maf0.1.rda")

#neighborDists <- neighborDists[keepIdxs]
#neighborDists <- neighborDists[-length(neighborDists)]
#neighborCors <- neighborCors[keepIdxs]
#neighborCors <- neighborCors[-length(neighborCors)]

