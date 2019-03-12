rm(list=ls())
library(data.table)
phenotype <- read.table("../data/phenotype/rawPhenotypeDatasets/Setaria_stomata_2016_GWAS_phenotype.csv",sep=",",header=TRUE,stringsAsFactors = FALSE)

#####Take mean of replicates in wet and dry####
widePheno <- data.table::dcast(setDT(phenotype),Genotype~Treatment,value.var=c("stoden_raw","stoden_pred"),fun.aggregate=mean,na.rm=T)
write.table(widePheno,"../data/phenotype/9.StomatalDensity.phenotypes.csv",sep=",",col.names=T,row.names=F)
