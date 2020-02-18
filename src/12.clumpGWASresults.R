#Expects intput to tell it which chromsome to clump#
args <- commandArgs(trailingOnly = TRUE)
chr <- ifelse(!(is.na(args[1])),args[1],4) #default to 4 for testing
message("Running clumping algorithm on chromosome ",chr)
####Make clumps of all vs all at n window size####
###Begins clumping at highest p-value SNPs first (so they are centered in a clump)
#####Go into each of those clumps and calculate LD between all SNPs in each clump####
######This should be similar, but not identical to PLINKs --clump command, in that I don't use 'index' SNPs to seed the clumps (except maybe for the stringent
##vs all)
#rm(list=ls())
source("ldCalculationFunction.R")
library(feather)
library(data.table)
library(plyr)
library(reshape2)
library(dplyr)
ldWin <- 250000 #distance around target SNP to consider SNPs in a clump, 250k is PLINK's default
r2 <- 0.75 #Threshold to include SNPs in a clump, 0.5 is PLINK's default

#load("../GWASresults/11.CombinedMLMMResultsWithFisherPval.filtered.rda") #loads object filterRes
#load("../GWASresults/11.allLTphenos.IB005to8.CombinedMLMMResultsWithFisherPval.filtered.rda") #loads object filterRes
snpTable <- feather::read_feather("../GWASresults/11.allLTphenos.IB003to4.CombinedMLMMResultsWithFisherPval.allSNPs.feather")
#use next to lines to further filter filterRes if desired
#fisherPcols <- grep("fisherP",colnames(filterRes),value=TRUE)
#filterRes <- filterRes[filterRes[, Reduce(`|`, lapply(.SD, `<=`, 1e-3)),.SDcols = fisherPcols]] #keep phenotype columns where fisherP <= firstThresh
#snpTable <- filterRes
#rm(filterRes)

snpTable$chrom <- as.numeric(sapply(strsplit(snpTable$SNP,"_"),"[",1))
snpTable$pos <- as.numeric(sapply(strsplit(snpTable$SNP,"_"),"[",2))

# snpTable <- fread("../GWASresults/11.GWASresults.0304_setaria_exp.ALL_days.BLUPS.csv",sep=",",stringsAsFactors = FALSE)
# snpTable2 <- fread("../GWASresults/11.GWASresults.050708_setaria_exp.ALL_days.BLUPS.csv",sep=",",stringsAsFactors = FALSE)
# snpTable$exp <- "0304"
# snpTable2$exp <- "050708"
# snpTable <- rbind(snpTable,snpTable2)
# rm(snpTable2)
# snpTable$Treatment <- sapply(strsplit(sapply(strsplit(snpTable$trait,"_"),"[",1),"\\."),"[",1)
# snpTable$DAP <- sapply(strsplit(sapply(strsplit(snpTable$trait,"_"),"[",1),"\\."),"[",2)
# snpTable$Phenotype <- sapply(strsplit(snpTable$trait,"\\d_"),"[",2)
# snpTable$chrom <- as.numeric(sapply(strsplit(snpTable$SNP,"_"),"[",1))
# snpTable$pos <- as.numeric(sapply(strsplit(snpTable$SNP,"_"),"[",2))

load("../data/genotype/5.filteredSNPs.noHighCorSNPs.2kbDistThresh.0.5neighborLD.0.975LDfilter.rda")#loads filterHighGeno, neighbors, filterHighInfo, filterHighResults

###Use phenotype to just get down to lines present in this trait###
#phenotype <- read.table("../data/phenotype/0304_setaria_exp.ALL_days.BLUPS.fromCharles.csv",sep=",",header=TRUE,stringsAsFactors = FALSE)
phenotype <- read.table("../data/phenotype/ALL_setaria_exp.ALL_days.TRANS_STAND.wet_dry_RATIO.csv",sep=",",header=TRUE,stringsAsFactors = FALSE)
colnames(phenotype)[1] <- "Genotype"
phenotypedGenos <- unique(phenotype$Genotype)

#phenotype <- read.table("../data/phenotype/050708_setaria_exp.ALL_days.BLUPS.fromCharles.csv",sep=",",header=TRUE,stringsAsFactors = FALSE)
#get to 1 column per treatment/trait
#meltPhenotype <- reshape2::melt(phenotype,id.vars=c("genotype","trt.day"))
#meltPhenotype <- meltPhenotype[which(!(is.na(meltPhenotype$value))),]
#phenotype <- reshape2::dcast(meltPhenotype,...~trt.day+variable)
#colnames(phenotype)[1] <- "Genotype"
#phenotypedGenos <- unique(c(phenotypedGenos,phenotype$Genotype))

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

row.names(filterHighGeno) <- paste(filterHighInfo$chr,filterHighInfo$pos,sep="_")
filterHighInfo$SNP <- paste(filterHighInfo$chr,filterHighInfo$pos,sep="_")

####Only keep SNPs that are in results###
genotype <- filterHighGeno[row.names(filterHighGeno) %in% snpTable$SNP,]
snpInfo <- filterHighInfo[filterHighInfo$SNP %in% snpTable$SNP,]
setnames(snpInfo,old = "chr", new = "chrom")
rm(filterHighGeno,filterHighInfo)
##Work flow
####Start at first SNP on the chromosome
#####Are there SNPs +/- dist 
######If yes, calculate LD between them
#######Remove any SNPs not in LD > r2 (removed SNPs go back into index SNP array)

clump <- 1 ####Current clump being analyzed
snpTable$clump <- NA
snpTable <- as.data.frame(snpTable)
snpTable <- snpTable[order(snpTable$avg.sv.height_wet_18_fisherP,as.numeric(snpTable$chrom),as.numeric(snpTable$pos)),]
#for(i in sort(unique(snpTable$chrom))){
message("Beginning clumping of chromosome ", chr)
thisChr <- snpTable[which(snpTable$chrom==chr),]
rm(snpTable)
thisGeno <- genotype[snpInfo$chrom==chr,]
thisInfo <- snpInfo[snpInfo$chrom==chr,]
for(j in 1:nrow(thisChr)){
  if(j %% 2500 == 0){
    cat("Analysing input row:",j,"of:",nrow(thisChr),"\n")
  }
  if(!(is.na(thisChr$clump[j]))){ #this snp is already in a clump
    next;
  }
  closeSnps <- thisChr[as.numeric(thisChr$pos) > as.numeric(thisChr$pos[j])-ldWin & 
                         as.numeric(thisChr$pos) < as.numeric(thisChr$pos[j])+ldWin,]
  closeSnps <- closeSnps[which(is.na(closeSnps$clump)),]
  if(nrow(closeSnps)==1){####Nothing around this SNP, it's a single snp clump####
    #snpTable$clump[snpTable$SNP==closeSnps$SNP] <- clump
    thisChr$clump[thisChr$SNP==closeSnps$SNP] <- clump
    clump <- clump+1      
  }else{
    ############Calculate LD between current SNP and others in window##############
    results <- data.frame()
    for(k in closeSnps$SNP){
      thisLD <- LDnumGeno(thisGeno[as.character(thisChr$SNP[j]),],thisGeno[k,])
      results <- rbind(results,data.frame(candidate=as.character(thisChr$SNP[j]),target=k,candidateBP=as.numeric(thisChr$pos[j]),dp=thisLD$`D'`,dpval=thisLD$`P-value`,r2=thisLD$`R^2`,n=thisLD$`n`,bp=closeSnps$pos[closeSnps$SNP==k][1]))
    }
    results <- results[results$r2 >= r2,]
    #snpTable$clump[snpTable$SNP %in% results$target] <- clump
    thisChr$clump[thisChr$SNP %in% results$target] <- clump
    clump <- clump+1
  }
}
#}

snpTable <- thisChr[order(as.numeric(thisChr$chrom),as.numeric(thisChr$pos)),]
rm(thisChr)
dir.create("../GWASresults/clumpRes", showWarnings = FALSE)
feather::write_feather(snpTable,path=paste0("../GWASresults/clumpRes/12.IB003toIB004.CombinedMLMMResultsWithFisherPval.allSNPs.clumpCol.chr",chr,".feather"))

####Notes from initial testing####
###firstThresh of 1e-3 returns final of 253 SNPs > logp5
###firstThresh of 1e-2 returns final of 360 SNPs > logp5
###firstThresh of 1e-3 returns final of 144 SNPs > logp6
###firstThresh of 1e-2 returns final of 186 SNPs > logp6
#The top 10 SNPs using both 1e-2 and 1e-3 filters are the same, although ranges vary slightly

#write.table(snpTable,"../results/12.clumpColumn.GWASresults.0304_setaria_exp.ALL_days.BLUPS.csv",sep=",",col.names=TRUE,row.names=FALSE)
#write.table(filterTable,"../results/12.1reppertrait.clumpColumn.GWASresults.0304_setaria_exp.ALL_days.BLUPS.csv",sep=",",col.names=TRUE,row.names=FALSE)

#This line gets just 'clumps' with two or more SNPs, don't think we want that yet
#snpTable[snpTable$clump %in% names(table(snpTable$clump)[table(snpTable$clump)>1]),]


# clumps <- snpTable[snpTable$clump %in% snpTable$clump[which(duplicated(snpTable$clump))],]
# clumps <- clumps[order(clumps$clump),]
# #####Break into clumps with two experiment types overlap#####
# typeOverlap <- ddply(clumps,.(clump),function(x){if(nrow(unique(as.data.frame(x[,"experiment"])))>=2){x}else{x[0,]}})
# yearOverlap <- ddply(typeOverlap,.(clump),function(x){if(nrow(unique(as.data.frame(x[,"year"])))>=2){x}else{x[0,]}})
# 
# #######Work out a scoring system for 'clumps'#########
# clumps$score <- 0
# for(i in unique(clumps$clump)){
#   thisClump <- clumps[clumps$clump==i,]
#   score <- 0
#   ###Are there any stringent SNPs###
#   #####If yes, check if SNPs in clump are from different years or from different types (e.g. root or ion)###
#   if(any(thisClump$model=="Stringent")){
#     if(length(unique(thisClump$year))>1 | length(unique(thisClump$type))>1){
#       score <- score + 12
#     }
#   }
#   ###Are there more than 2 PC traits###
#   pcThisClump <- thisClump[grep("PC",thisClump$phenotype),]
#   # if(nrow(pcThisClump)>1){
#   #   if(length(unique(pcThisClump$year))>1 | length(unique(pcThisClump$type))>1){
#   #     score <- score + 9
#   #   }
#   # }else 
#   if(nrow(pcThisClump)>1){
#     if(length(unique(thisClump$year))>1 | length(unique(thisClump$type))>1){
#       score <- score + 8
#     }
#   }
#   ##Same trait, different years
#   if(any(rowSums(table(thisClump$phenotype,thisClump$year))>1)){
#     if(any(thisClump$model=="stringent")){
#       score <- score + 10
#     }else{
#       score <- score + 7
#     }
#   }
#   if(length(unique(thisClump$type))>1){
#     if(any(thisClump$model=="stringent")){
#       score <- score + 9
#     }else{
#       score <- score + 7
#     }
#   }
#   clumps$score[clumps$clump==i] <- score
# }
# 
# clumpTable <- data.frame()
# for(i in unique(clumps$clump)){
#   thisClump <- clumps[clumps$clump==i,]
#   pcThisClump <- thisClump[grep("PC",thisClump$phenotype),]
#   clumpTable <- rbind(clumpTable,data.frame(
#     clump = i,
#     NumSNPs = nrow(thisClump),
#     NumStringentSNPs = length(which(thisClump$model=="Stringent")),
#     NumPCtraits = length(grep("PC",thisClump$phenotype)),
#     StringentSNPOverlapsAnotherType = ifelse(any(thisClump$model=="Stringent") & length(unique(thisClump$type))>1,TRUE,FALSE),
#     StringentSNPOverlapsAnotherYear = ifelse(any(thisClump$model=="Stringent") & length(unique(thisClump$year))>1,TRUE,FALSE),
#     MultiplePCoverlapType = ifelse(length(unique(pcThisClump$type))>1,TRUE,FALSE),
#     MultiplePCoverlapYear = ifelse(length(unique(pcThisClump$year))>1,TRUE,FALSE),
#     SameTraitDifferentYears = ifelse(any(rowSums(table(thisClump$phenotype,thisClump$year))>1),TRUE,FALSE),
#     DifferentTypes = ifelse(length(unique(thisClump$type))>1,TRUE,FALSE),
#     GelAndFieldTrait = ifelse(any(thisClump$year=="Gel") & any(thisClump$year != "Gel"),TRUE,FALSE)
#   ))
# }
# 
# write.table(clumps,"../results/10.results.clumps.csv",sep=",",row.names=FALSE,col.names=TRUE)
# write.table(clumpTable,"../results/10.results.clumpInfotable.csv",sep=",",row.names=FALSE,col.names = TRUE)
# #write.table(clumps[clumps$clump %in% c(1723,1775,2641,2801,1691,634),],"../results/10.intersetingClumps.csv",sep=",",row.names=FALSE,col.names=TRUE)
# #10 pts - stringent SNP with a SNP from another year or type (e.g. ion with root)
# # 9 pts - PC SNPs found across years or types (e.g. ion with root)
# # 8 pts - PC SNP found with a nonPC SNP across years or types
# # 7 pts - Same SNP found for same trait in multiple years
# # 6 pts - SNP of different 'types' found
# ###Gel overlaps with?###
# # 
# # rootSnps <- olClumps[olClumps$experiment %in% c("dirt","gel","xray"),]
# # rootSnps <- rootSnps[rootSnps$clump %in% rootSnps$clump[which(duplicated(rootSnps$clump))],]
# # rootSnpsOverlapDiffYears <- ddply(rootSnps,.(clump),function(x){if(nrow(unique(as.data.frame(x[,"year"])))>=2){x}else{x[0,]}})
# # 
# # #####overlaps for ion across years####
# # ionSnps <- olClumps[olClumps$experiment=="ion",]
# # ionSnps <- ionSnps[ionSnps$clump %in% ionSnps$clump[which(duplicated(ionSnps$clump))],]
# # ionSnpsOverlapDiffYears <- ddply(ionSnps,.(clump),function(x){if(nrow(unique(as.data.frame(x[,"year"])))>=2){x}else{x[0,]}})
# # 
# ######Point system to code traits?
# ###Across years=x points, across phenotypes=y points...
# #####Of loci in best is there any evidence in others....
# #####Gel and then gel+support in field experiment....
# #####Then multi year root trait for any datapoint
# #############
# ###Summary table of loci, is it in gel and a field
# ####If yes, is it in selective or all
# ###What is the LD, how fast does it decay, what are the genes?
# # write.table(olIonRootOverlap,"../results/11.LDbasedIonRootOverlap.csv",sep=",",row.names=FALSE,col.names=TRUE)
# # write.table(rootSnpsOverlapDiffYears,"../results/11.LDbasedRootOverlapDiffYears.csv",sep=",",row.names=FALSE,col.names=TRUE)
# # write.table(ionSnpsOverlapDiffYears,"../results/11.LDbasedIonOverlapDiffYears.csv",sep=",",row.names=FALSE,col.names=TRUE)
# # write.table(olClumps,"../results/11.LDbasedAllOverlaps.csv",sep=",",row.names=FALSE,col.names=TRUE)
# # 
# # 
