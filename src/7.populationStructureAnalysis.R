rm(list)
library(data.table)
library(ggplot2)
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

####Look at kinship of lemnatech lines for very similar lines

