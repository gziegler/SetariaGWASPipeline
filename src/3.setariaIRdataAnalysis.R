rm(list=ls())
library(lme4)
library(ggplot2)
library(data.table)
library(ggpubr)
library(lmomco)
phenotype <- read.table("../data/phenotype/rawPhenotypeDatasets/Setaria_IR_2016_datsetset_GWAS.csv",header=TRUE,stringsAsFactors = FALSE,sep=",")
table(phenotype$Awning[phenotype$Genotype=="B100"])
#ggplot(phenotype,aes(y=Jul_11_IR,x=factor(Awning),fill=Treatment))+geom_boxplot()
#ggplot(phenotype,aes(y=Jul_21_IR,x=factor(Awning),fill=Treatment))+geom_boxplot()
ggplot(phenotype[phenotype$Genotype=="B100",],aes(y=Jul_11_IR,x=factor(Awning),fill=Treatment))+geom_boxplot()+theme_bw()+theme(axis.text.x = element_text(size=22),axis.text.y=element_text(size=22))
ggplot(phenotype[phenotype$Genotype=="B100",],aes(y=Jul_16_IR,x=factor(Awning),fill=Treatment))+geom_boxplot()+theme_bw()+theme(axis.text.x = element_text(size=22),axis.text.y=element_text(size=22))
ggplot(phenotype[phenotype$Genotype=="B100",],aes(y=Jul_21_IR,x=factor(Awning),fill=Treatment))+geom_boxplot()+theme_bw()+theme(axis.text.x = element_text(size=22),axis.text.y=element_text(size=22))

widePheno <- data.table::dcast(setDT(phenotype),Genotype~Treatment,value.var=c("Jul_11_IR","Jul_16_IR","Jul_21_IR","Awning"),fun.aggregate=mean,na.rm=T)

phenotype$Check <- ifelse(phenotype$Genotype=="B100",0,1)
phenotype$GenoCheckOrTest <- ifelse(phenotype$Genotype=="B100",phenotype$Genotype,999)
phenotype[,c("Check","GenoCheckOrTest")]
#lmer(Jul_16_IR~Check)
#Appears to be an awning effect

#ggplot(widePheno,aes(x = Jul_11_IR_dry, y=Jul_21_IR_dry))+geom_point()

traits <- colnames(phenotype)[6:8]
H2 <- data.frame()
for(i in traits){
  thisTrait <- phenotype[,c("Genotype","Treatment","Awning","Check","GenoCheckOrTest",i),with=FALSE]
  thisTrait <- thisTrait[complete.cases(thisTrait),]
  thisTrait$Awning <- factor(thisTrait$Awning)
  thisTrait$Treatment <- factor(thisTrait$Treatment)
  thisTrait$Genotype <- factor(thisTrait$Genotype)
  thisTrait$GenoCheckOrTest <- factor(thisTrait$GenoCheckOrTest)
  thisTrait$Check <- factor(thisTrait$Check)
  #####This is for an attempt at analysis as an 'augmented block design'
  ###but I think it's more for testing if there is significant differences between test plots vs check plots
  #####so i didn't end up using it. See https://passel.unl.edu/pages/index.php?alllessons=1
  ##and https://articles.extension.org/pages/60430/introduction-to-the-augmented-experimental-design-webinar
  #augModel <- lmer(Jul_16_IR ~ 0 + Treatment + GenoCheckOrTest + (1|Awning) + (1|Genotype:Check),data=thisTrait)
  
  thisForm0 <- formula(paste(i,"~ 0 + (1|Genotype)"))
  thisForm1 <- formula(paste(i,"~ 0 + Treatment + (1|Genotype)"))
  thisForm2 <- formula(paste(i,"~ 0 + Treatment + (1|Awning) + (1|Genotype)"))
  thisForm3 <- formula(paste(i,"~ 0 + Treatment + (1|Awning) + (1|Genotype) + (1|Genotype:Treatment)"))
  #m4 is the model Jeff recommends, should be very similar to m2 except m2 will shrink Genotype to the mean in cases of uncertainty?
  thisForm4 <- formula(paste(i,"~ 0 + (1|Awning) + Genotype:Treatment"))
  m0 <- lmer(thisForm0,data=thisTrait)
  m1 <- lmer(thisForm1,data=thisTrait)
  m2 <- lmer(thisForm2,data=thisTrait)
  m3 <- lmer(thisForm3,data=thisTrait)
  m4 <- lmer(thisForm4,data=thisTrait)
  print(anova(m0,m1,m2,m3,m4))
  
  m4genoBLUEs <- as.data.frame(fixef(m4))
  colnames(m4genoBLUEs) <- paste0(i,"_BLUE")
  
  m4genoBLUEs$Genotype <- sub("Genotype","",sapply(strsplit(row.names(m4genoBLUEs),":"),"[",1))
  m4genoBLUEs$Treatment <- sub("Treatment","",sapply(strsplit(row.names(m4genoBLUEs),":"),"[",2))
  widem4 <- data.table::dcast(setDT(m4genoBLUEs),Genotype~Treatment,value.var=paste0(i,"_BLUE"),fun.aggregate=mean,na.rm=T)
  colnames(widem4)[2:3] <- paste(i,colnames(widem4)[2:3],"BLUE",sep="_")
  widePheno <- merge(widePheno,widem4,all=T)
  
  ran2 <- ranef(m2)
  m2genoBLUPs <- ran2$Genotype
  m2genoBLUPs <- data.frame(Genotype=row.names(m2genoBLUPs),BLUP=m2genoBLUPs$`(Intercept)`)
  colnames(m2genoBLUPs)[2] <- paste0(i,"_BLUP")
  widePheno <- merge(widePheno,m2genoBLUPs,all=T)
  
  ####Calculate heritability from m2, but without B100 checks
  thisTraitNoCheck <- thisTrait[thisTrait$Genotype != "B100",]
  ###Remove non-
  m2 <- lmer(thisForm2,data=thisTraitNoCheck)
  #first get varaince of fixed treatment effect
  ###There seem to be many ideas of how to get this, each gives similar but not identical variance
  #####1. is to fit model with and without Treatment and subtract the variance
  ########this is the method brianna used
  #####2. fit a linear model with just treatment and get the variance of the predicted values
  #####3. Get variance of fixed effects by multiplying coefficients by design matrix m2@pp$X is same as model.matrix(m2)
  ####See https://stats.stackexchange.com/questions/136087/lmm-how-do-i-calculate-a-standard-deviation-on-the-variance-explained-by-fixed?rq=1
  ####and https://jonlefcheck.net/2013/03/13/r2-for-linear-mixed-effects-models/ and 
  ####https://github.com/jslefche/piecewiseSEM/blob/master/R/rsquared.R
  ##more discussion here:http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#how-do-i-compute-a-coefficient-of-determination-r2-or-an-analogue-for-glmms
  #####4. Fit a model where everything is random effects, this give the highest variance for Treatment
  ###Fixed Ef Variance using step 1 option
  m2woTreat <- update(m2,.~.-Treatment)
  re <- as.data.frame(VarCorr(m2))$vcov
  
  reWT <- as.data.frame(VarCorr(m2woTreat))$vcov
  fixefVarDiffBetweenModels <- sum(reWT)-sum(re)
  fixefVarLM <- var(predict(lm(formula(paste0(i,"~Treatment")),data = thisTraitNoCheck))) 
  fixefVarModelMatrix <- var(as.vector(lme4::fixef(m2) %*% t(m2@pp$X)))
  fixefVarRandModel <- as.data.frame(VarCorr(update(m2,.~.-Treatment+(1|Treatment))))$vcov[3]
  
  ###The fixefVarModelMatrix makes the most sense to me and it's essentially giving the same thing as fixefVarLM, so I'm going with that answer
  ####I'm not sure why the fixefVarRandModel option give such a large variance 
  ######perhaps because Awning and Treatment are correlated (also sum of variance is much larger than var(phenotype$Jul_21_IR))
  tot.var <- sum(fixefVarModelMatrix,re)
  
  #Get reps in each condition to weight the variance by because we don't have fully replicated design (see convo with Lipka/Feldman)
  reps.t1<-as.character(unique(thisTraitNoCheck$Genotype[thisTraitNoCheck$Treatment == "dry"]))
  reps.t2<-as.character(unique(thisTraitNoCheck$Genotype[thisTraitNoCheck$Treatment == "wet"]))

  unique.combined <- c(as.character(reps.t1), as.character(reps.t2))
  
  freq.unique.combined <- table(unique.combined)
  
  # Calculate the harmonic mean replication within treatment blocks
  hm_treatment<-harmonic.mean(freq.unique.combined)$harmean
  
  reps.total<-table(thisTraitNoCheck$Genotype)
  hm_total<-harmonic.mean(reps.total)$harmean
  
  betweenTreatmentH2 <- re[1]/(re[1]+re[2]+fixefVarModelMatrix+re[3]/hm_total)
  H2 <- rbind(H2,data.frame(Phenotype=i,Treatment="Both",H2=betweenTreatmentH2,nLines=length(unique(phenotype$Genotype[which(!(is.na(phenotype[,..i])))]))))
  #########Calculate within treatment H2#####
  for(j in c("wet","dry")){
    withinForm <- formula(paste(i,"~ 1 + (1|Awning) + (1|Genotype)"))
    thisPheno <- thisTraitNoCheck[thisTraitNoCheck$Treatment==j,]
    thisPheno <- droplevels(thisPheno)
    withinMod <- lmer(withinForm,data=thisPheno)
    re <- as.data.frame(VarCorr(withinMod))$vcov
    withinTreatmentH2 <- re[1]/(sum(re))
    H2 <- rbind(H2,data.frame(Phenotype=i,Treatment=j,H2=withinTreatmentH2,nLines=length(unique(phenotype$Genotype[which(!(is.na(phenotype[,..i])))]))))
  }
  
}

#Due to most lines not having sufficient replicates, H2 results aren't valid for Jul_16_IR
H2 <- H2[H2$Phenotype != "Jul_16_IR",]
#H2$Phenotype <- factor(H2$Phenotype,levels=c("Jul_11_IR","Jul_16_IR","Jul_21_IR"))
H2$Phenotype <- factor(H2$Phenotype,levels=c("Jul_11_IR","Jul_21_IR"))
ggplot(H2,aes(x=Phenotype,fill=Treatment,y=H2))+geom_col(position="dodge")+theme_bw()+theme(axis.text.x = element_text(size=22),axis.text.y=element_text(size=22))

####Compare actual values to values from m4 fixed effect (BLUEs)
widePheno$Awning_dry <- as.factor(widePheno$Awning_dry)
widePheno$Awning_wet <- as.factor(widePheno$Awning_wet)
pdf("../results/2.setariaIR.BLUEsAndBlupsVsOriginal.pdf",width=12,height=8)
for(i in traits){
  for(j in c("wet","dry")){
    sp <- ggscatter(widePheno, x = paste(i,j,sep="_"), y = paste(i,j,"BLUE",sep="_"),
                    color = paste("Awning",j,sep="_"),
                    add = "reg.line",  # Add regression line
                    title= i,
                    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE # Add confidence interval
    )
    # Add correlation coefficient
    sp <- sp + stat_cor(method = "pearson")
    
    print(sp)
  }
}

for(i in traits){
  for(j in c("wet","dry")){
    sp <- ggscatter(widePheno, x = paste(i,j,sep="_"), y = paste(i,"BLUP",sep="_"),
                    color = paste("Awning",j,sep="_"),
                    add = "reg.line",  # Add regression line
                    title= i,
                    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE # Add confidence interval
    )
    # Add correlation coefficient
    sp <- sp + stat_cor(method = "pearson")
    
    print(sp)
  }
}

for(i in traits){
  for(j in c("wet","dry")){
    sp <- ggscatter(widePheno, x = paste(i,j,"BLUE",sep="_"), y = paste(i,"BLUP",sep="_"),
                    color = paste("Awning",j,sep="_"),
                    add = "reg.line",  # Add regression line
                    title= i,
                    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE # Add confidence interval
    )
    # Add correlation coefficient
    sp <- sp + stat_cor(method = "pearson")
    
    print(sp)
  }
}

for(i in traits){
  sp <- ggscatter(widePheno, x = paste(i,"dry","BLUE",sep="_"), y = paste(i,"wet","BLUE",sep="_"),
                  color = paste("Awning",j,sep="_"),
                  add = "reg.line",  # Add regression line
                  title= i,
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE # Add confidence interval
  )
  # Add correlation coefficient
  sp <- sp + stat_cor(method = "pearson")
  
  print(sp)
  
  sp <- ggscatter(widePheno, x = paste(i,"dry",sep="_"), y = paste(i,"wet",sep="_"),
                  color = paste("Awning",j,sep="_"),
                  add = "reg.line",  # Add regression line
                  title= i,
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE # Add confidence interval
  )
  # Add correlation coefficient
  sp <- sp + stat_cor(method = "pearson")
  
  print(sp)
  
}

dev.off()

sp <- ggscatter(widePheno, x = "Jul_11_IR_wet_BLUE", y = "Jul_16_IR_wet_BLUE",
                #color = paste("Awning",j,sep="_"),
                add = "reg.line",  # Add regression line
                title= i,
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp <- sp + stat_cor(method = "pearson")
print(sp)

sp <- ggscatter(widePheno, x = "Jul_11_IR_BLUP", y = "Jul_16_IR_BLUP",
                #color = paste("Awning",j,sep="_"),
                add = "reg.line",  # Add regression line
                title= i,
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp <- sp + stat_cor(method = "pearson")
print(sp)

write.table(widePheno,"../data/phenotype/2.Setaria_IR_2016_datsetset_GWAS.BLUPsandBLUEs.csv",sep=",",row.names=FALSE,col.names=T)
write.table(H2,"../results/2.setariaIR.H2.csv",sep=",",row.names = FALSE,col.names=T)

#######Open csv containing info about lines in Genotype file#####
genoInfo <- read.table("../data/genotype/Setaria_597_diversity_samples.csv",sep=",",header=TRUE,stringsAsFactors = FALSE,comment.char = "")
genoInfo$Genotype <- gsub("_setaria_12","",genoInfo$New_name)
length(intersect(genoInfo$Genotype,widePheno$Genotype))
setdiff(widePheno$Genotype,genoInfo$Genotype)
####Just 2 (3 including b100) genotypes not found in genotype file
##writing out table of names to keep for vcf filtering
keepLines <- as.data.frame(sapply(genoInfo$LIB[genoInfo$Genotype %in% c("A10.1",intersect(genoInfo$Genotype,widePheno$Genotype))],function(x){paste(rep(x,4),collapse = "_")}))
write.table(keepLines,"../data/genotype/keepLines.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)

###m2 seems to fit the model the best and indicates not a lot of GxE in this population
####Fitting Genotype:Treatment as a fixed effect may be easier to interpret if not interested in 'Global' BLUPs
#####This also prevents the 'shrinkage' of Genotypic values
