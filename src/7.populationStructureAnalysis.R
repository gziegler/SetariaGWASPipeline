rm(list)
library(data.table)
library(ggplot2)
library(gplots)
load("../data/genotype/6.AstleBalding.synbreed.kinship.rda")
load("../data/genotype/6.Eigenstrat.population.structure.50PCs.rda")

#######Open csv containing info about lines in Genotype file#####
genoInfo <- read.table("../data/genotype/Setaria_597_diversity_samples.csv",sep=",",header=TRUE,stringsAsFactors = FALSE,comment.char = "")
genoInfo$Genotype <- gsub("_setaria_12","",genoInfo$New_name)

#ID is repped 4 times in geno files
genoInfo$from <- sapply(genoInfo$LIB,function(x){paste(rep(x,4),collapse = "_")})

commonLines <- intersect(genoInfo$from,colnames(kinship))
genoInfo <- genoInfo[genoInfo$from %in% commonLines,]
kinship <- kinship[commonLines,commonLines]
structData <- structData[commonLines,]

genoInfo <- genoInfo[order(match(genoInfo$from,colnames(kinship))),]
stopifnot(identical(genoInfo$from,colnames(kinship)))
stopifnot(identical(genoInfo$from,row.names(structData)))

colnames(kinship) <- genoInfo$Genotype
row.names(kinship) <- genoInfo$Genotype
row.names(structData) <- genoInfo$Genotype

####Read in lines used on phenotype
lemnaLines <- read.table("../data/genotype/setaria_bellweather_2017.csv",sep=",",header=TRUE,stringsAsFactors = FALSE)
#PH_023 is duplicated, so there are only 293 lemnatech lines
lemnaLines[duplicated(lemnaLines$ID),]

genoInfo$lemnaID <- NA
genoInfo$lemnaID[genoInfo$New_name %in% intersect(genoInfo$New_name,lemnaLines$ID)] <- genoInfo$New_name[genoInfo$New_name %in% intersect(genoInfo$New_name,lemnaLines$ID)]

genoInfo$lemnaID[grepl("PH",genoInfo$Name) & is.na(genoInfo$lemnaID)] <- gsub("_2017","",gsub("PH_2014_|PH_2016_","PH_",genoInfo$Name[grepl("PH",genoInfo$Name) & is.na(genoInfo$lemnaID)]))
genoInfo$lemnaID[!(genoInfo$lemnaID %in% intersect(genoInfo$lemnaID,lemnaLines$ID))] <- NA

genoInfo$lemnaID[grep("^MF\\d+_2017",genoInfo$Name)] <- gsub("_2017","",gsub("MF","MF_",genoInfo$Name[grep("^MF\\d+_2017",genoInfo$Name)]))

#######Can't match 1 line: TB_setaria_12_0048
setdiff(lemnaLines$ID,genoInfo$lemnaID)

LTlines <- genoInfo[genoInfo$lemnaID %in% intersect(lemnaLines$ID,genoInfo$lemnaID),]

########Make structure plot colored by lines run on phenotyper
#add a column to structure file with if in lemnatech or not
structTable <- as.data.frame(structData)
colnames(structTable) <- gsub("V","PC",colnames(structTable))
structTable$Genotype <- row.names(structTable)
structTable$InLemnaTec <- "No"

structTable$InLemnaTec[structTable$Genotype %in% LTlines$Genotype] <- "Yes"

LTlines <- merge(LTlines,lemnaLines[,c("ID","percGerm")],by.x="lemnaID",by.y="ID")

structTable$InLemnaTec <- "No"
structTable$InLemnaTec[structTable$Genotype %in% LTlines$Genotype] <- "Yes"

structTable <- merge(structTable,LTlines[,c("Genotype","percGerm","lemnaID")],by="Genotype",all.x=T)
structTable$percGermOrig <- structTable$percGerm
structTable$percGerm[structTable$percGerm<=0.5] <- "<=0.5"
structTable$percGerm[which(as.numeric(structTable$percGerm)>0.5)] <- ">0.5"

structTable$percGerm[which(is.na(structTable$percGerm))] <- "NP"

####Read in pop structure from Pu
Qmat <- read.table("../data/genotype/Qmatrix.fastSTRUCTURE.fromPU.tab",header=FALSE,stringsAsFactors = FALSE)
colnames(Qmat) <- c("Genotype","Pop1","Pop2","Pop3","Pop4")
Qmat$Pop <- apply(Qmat,1,function(x){paste0("Pop",which.max(x[2:5]))})
Qmat$Pop <- apply(Qmat,1,function(x){ifelse(any(x[2:5]>0.6),x[6],"Admixed")})

structTable <- merge(structTable,Qmat,by.x="lemnaID",by.y="Genotype",all.x=T,all.y=F)

pdf("../results/7.popStructure.totalVs2017lemnateclines.pdf",width=12,height=9)
print(ggplot(structTable,aes(x=PC1,y=PC2,color=InLemnaTec))+geom_point()+theme_bw())
print(ggplot(structTable,aes(x=PC1,y=PC2,color=percGerm))+geom_point()+theme_bw())
print(ggplot(structTable,aes(x=PC1,y=PC2,color=Pop))+geom_point()+theme_bw())
dev.off()

write.table(structTable,"../data/genotype/7.popStructureMergedwithLemnatecInfo.csv",sep=",",row.names=FALSE,col.names=T)

####Look at kinship of lemnatech lines for very similar lines
kinDF <- as.data.frame(kinship)
kinDF$Genotype <- row.names(kinDF)
kinDF <- kinDF[kinDF$Genotype %in% intersect(kinDF$Genotype,LTlines$Genotype),c("Genotype",intersect(kinDF$Genotype,LTlines$Genotype))]
LTlines <- LTlines[which(!(duplicated(LTlines$Genotype))),]
meltKin <- melt(kinDF)
meltKin <- merge(meltKin,LTlines[,c("Genotype","percGerm","lemnaID")],by="Genotype",all.x=T)
meltKin <- merge(meltKin,LTlines[,c("Genotype","percGerm","lemnaID")],by.x="variable",by.y="Genotype",all.x=T)
meltKin$germDiff <- meltKin$percGerm.x-meltKin$percGerm.y
meltKin <- meltKin[meltKin$germDiff>0,]
meltKin <- meltKin[meltKin$value>=2,]
meltKin <- meltKin[order(meltKin$germDiff,decreasing = T),]
length(unique(meltKin$lemnaID.y[meltKin$germDiff>0.5]))
hist(meltKin$value[meltKin$germDiff>0.5])
discardKin <- meltKin[meltKin$germDiff>0.5,]
#discardKin <- discardKin[discardKin$percGerm.y<=.40,]
########List of lines more similar to each other than >90% of other lines, 
########but with a >50% difference in germination percentage to closely related lines
#######Ended up using the distance matrix to filter instead, but they both actually give pretty similar results

discardLines <- unique(discardKin$lemnaID.y)
structTable$discardLines <- structTable$InLemnaTec
structTable$discardLines[structTable$lemnaID %in% discardLines] <- "Discard"
structTable[which(structTable$discardLines=="Discard"),c(1:4,53:61)]
structTable$Label <- ""
structTable$Label[structTable$lemnaID %in% discardLines] <- structTable$Genotype[structTable$lemnaID %in% discardLines]
print(ggplot(structTable,aes(x=PC1,y=PC2,color=discardLines))+geom_point(size=3,alpha=0.5)+geom_text(aes(label=Label))+theme_bw())
discardKin[discardKin$variable=="TB_0333",]

kinship <- kinship[which(row.names(kinship) %in% LTlines$Genotype),which(colnames(kinship) %in% LTlines$Genotype)]
pdf(paste("../results/7.kinshipOfLemnatecLines.pdf",sep=""), width = 12, height = 12)
par(mar = c(25,25,25,25))
heatmap.2(kinship,  cexRow =.2, cexCol = 0.2, col=rev(heat.colors(256)), scale="none", symkey=FALSE, trace="none")
dev.off()

##try again but with distance matrix
distKin <- as.matrix(dist(kinship))
distDF <- as.data.frame(distKin)
distDF$Genotype <- row.names(distDF)
meltDist <- melt(distDF)
meltDist <- merge(meltDist,LTlines[,c("Genotype","percGerm","lemnaID")],by="Genotype",all.x=T)
meltDist <- merge(meltDist,LTlines[,c("Genotype","percGerm","lemnaID")],by.x="variable",by.y="Genotype",all.x=T)
meltDist$germDiff <- meltDist$percGerm.x-meltDist$percGerm.y
meltDist <- meltDist[meltDist$germDiff>0,]
meltDist <- meltDist[meltDist$value<2.0,]
meltDist <- meltDist[order(meltDist$germDiff,decreasing = T),]
diffCut <- 0.25
length(unique(meltDist$lemnaID.y[meltDist$germDiff>diffCut]))
hist(meltDist$value[meltDist$germDiff>diffCut])
discardDist <- meltDist[meltDist$germDiff>diffCut,]
discardLinesDist <- unique(discardDist$lemnaID.y)

structTable$discardLines <- structTable$InLemnaTec
structTable$discardLines[structTable$lemnaID %in% discardLinesDist] <- "Discard"
structTable[which(structTable$discardLines=="Discard"),c(1:4,53:61)]
structTable$Label <- ""
structTable$Label[structTable$lemnaID %in% discardLinesDist] <- structTable$Genotype[structTable$lemnaID %in% discardLinesDist]

pdf("../results/7.PC1vsPC2.discardedLinesLabeled.pdf",width=12,height=9)
print(ggplot(structTable,aes(x=PC1,y=PC2,color=discardLines))+geom_point(size=3,alpha=0.5)+geom_text(aes(label=Label))+theme_bw())
dev.off()

keepLines <- LTlines[which(!(LTlines$lemnaID %in% discardLinesDist)),]
write.table(keepLines,"../data/genotype/7.selectedLinesFor2018LemnaTec.csv",sep=",",col.names=T,row.names=F)
write.table(structTable,"../data/genotype/7.popStructureMergedwithLemnatecInfo.csv",sep=",",row.names=FALSE,col.names=T)

#discardKin[discardKin$variable=="TB_0333",]
#discardDist[discardDist$variable=="TB_0060",]

#structTable$Label[structTable$Genotype=="TB_setaria_16_0671"] <- "TB_setaria_16_0671"

# ####Some code for filtering but just with distance
# ##try again but with distance matrix
# distKin <- as.matrix(dist(kinship))
# distDF <- as.data.frame(distKin)
# distDF$Genotype <- row.names(distDF)
# meltDist <- melt(distDF)
# meltDist <- merge(meltDist,LTlines[,c("Genotype","percGerm","lemnaID")],by="Genotype",all.x=T)
# meltDist <- merge(meltDist,LTlines[,c("Genotype","percGerm","lemnaID")],by.x="variable",by.y="Genotype",all.x=T)
# meltDist$germDiff <- meltDist$percGerm.x-meltDist$percGerm.y
# meltDist <- meltDist[meltDist$germDiff>0,]
# meltDist <- meltDist[meltDist$value<1,]
# meltDist <- meltDist[meltDist$Genotype != meltDist$variable,]
# length(unique(meltDist$lemnaID.y))
# distRemove <- unique(meltDist$lemnaID.y)
# length(unique(c(meltDist$lemnaID.y,discardLinesDist)))
# 
# 
# structTable$discardLines[structTable$lemnaID %in% distRemove] <- "Discard by Similarity"
# structTable$Label[structTable$lemnaID %in% distRemove] <- structTable$Genotype[structTable$lemnaID %in% distRemove]
# print(ggplot(structTable,aes(x=PC1,y=PC2,color=discardLines))+geom_point(size=3,alpha=0.5)+geom_text(aes(label=Label))+theme_bw())
# structTable[structTable$discardLines %in% c("Discard","Discard by Similarity"),c(1:4,51:61)]
