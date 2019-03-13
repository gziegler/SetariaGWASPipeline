irRes <- read.table("../GWASresults/11.allStomataresults.csv",stringsAsFactors = FALSE,sep=",",header=TRUE)
irRes$chr <- sapply(strsplit(irRes$SNP,"_"),"[",1)
irRes$bp <- sapply(strsplit(irRes$SNP,"_"),"[",2)
write.table(irRes,"../GWASresults/12.allStomataresults.ZBrowseReady.csv",sep=",",col.names = TRUE,row.names=FALSE)
