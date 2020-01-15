####Heavily adapted from LD function in genetics library#######
#####Expects 0, 1, 2 genotypes, where 0 is the major allele#####
#####Load by sourcing this function####
#####LDnumGeno(g1vec,g2vec)#####
#Test examples
#g1num <- c(1, NA, 0, NA, 1, NA, 0, 1, 0, 0, 1, 2, 0, 1, 1, 0, NA, 1, 1, NA)
#g2num <- c(1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 2, 1, 2, 1, 0, 1, 1, 1, 2)
#g2num <- c(.5,.5,0,.5,0,.5,.5,.5,.5,0,.5,0,.5,.5,.5,0,.5,.5,.5,0)
#g2num01 <- c(.5, .5, 0, .5, 0, .5, .5, .5,.5, 0, .5, 1, .5, 1, .5, 0,.5, .5, .5, 1)
#g1num01 <- c(.5,    NA, 0,  NA, .5,    NA, 0, .5,0, 0, .5, 1, 0, .5, .5, 0,NA, .5, .5,   NA)

LDnumGeno <- function(g1num,g2num){
  calcAlleleFreq <- function(alleleVec){
    tot1 <- length(alleleVec[which(alleleVec==2)])*2
    tot0 <- length(alleleVec[which(alleleVec==0)])*2
    totHet <- length(alleleVec[which(alleleVec==1)])
    totAllele <- length(alleleVec[which(!(is.na(alleleVec)))])*2
    prop.A.num <- c((tot0+totHet)/totAllele,(tot1+totHet)/totAllele)
    names(prop.A.num) <- c(0,2)
    prop.A.num
  }
  prop.A <- calcAlleleFreq(g1num)
  prop.B <- calcAlleleFreq(g2num)
  
  major.A <- names(prop.A)[which.max(prop.A)]
  minor.A <- names(prop.A)[which.min(prop.A)]
  major.B <- names(prop.B)[which.max(prop.B)]
  minor.B <- names(prop.B)[which.min(prop.B)]

  pA <- max(prop.A, na.rm=TRUE)
  pB <- max(prop.B, na.rm=TRUE)
  pa <- 1-pA
  pb <- 1-pB

  Dmin <- max(-pA*pB, -pa*pb)
  pmin <- pA*pB + Dmin;

  Dmax <- min(pA*pb, pB*pa);
  pmax <- pA*pB + Dmax;

  #convert allele vector to 0,1,2,NA (0 is homo minor, 1 is het, 2 is homo major), then make a frequency table
  #  comparing the two alleles
  # g1num01[which(g1num01==major.A01)] <- 2
  # g1num01[which(g1num01==minor.A01)] <- 0
  # g1num01[which(g1num01==0.5)] <- 1
  # 
  # g2num01[which(g2num01==major.B01)] <- 2
  # g2num01[which(g2num01==minor.B01)] <- 0
  # g2num01[which(g2num01==0.5)] <- 1
  
  counts <- table(g1num,g2num)
  #n3x3 <- counts
  #make highest frequency values in upper left (this should mean 0's in upper left), which is what happens
  #if coded 0 is homo major...
  # n3x3 <- matrix(0, nrow=3, ncol=3)
  # colnames(n3x3) <- rownames(n3x3) <- 0:2
  # for(i in rownames(counts))
  #   for(j in colnames(counts))
  #     n3x3[3-as.numeric(i),3-as.numeric(j)] <- counts[i,j]
  n3x3 <- matrix(0, nrow=3, ncol=3)
  colnames(n3x3) <- rownames(n3x3) <- 0:2
  for(i in rownames(counts))
    for(j in colnames(counts))
      n3x3[as.numeric(i)+1,as.numeric(j)+1] <- counts[i,j]
  
  loglik <- function(pAB,...)
  {
    (2*n3x3[1,1]+n3x3[1,2]+n3x3[2,1])*log(pAB) +
      (2*n3x3[1,3]+n3x3[1,2]+n3x3[2,3])*log(pA-pAB) +
      (2*n3x3[3,1]+n3x3[2,1]+n3x3[3,2])*log(pB-pAB) +
      (2*n3x3[3,3]+n3x3[3,2]+n3x3[2,3])*log(1-pA-pB+pAB) + 
      n3x3[2,2]*log(pAB*(1-pA-pB+pAB) + (pA-pAB)*(pB-pAB))
  }
  
  # SAS code uses:
  #
  #s <- seq(pmin+0.0001,pmax-0.0001,by=0.0001)
  #lldmx <- loglik(s)
  #maxi <- which.max(lldmx)
  #pAB <- s[maxi]
  
  # but this should be faster:
  solution <- optimize(
    loglik,
    lower=pmin+.Machine$double.eps,
    upper=pmax-.Machine$double.eps,
    maximum=TRUE
  )
  pAB <- solution$maximum
  
  estD <- pAB - pA*pB
  if (estD>0)  
    {estDp <- estD / Dmax
  }else{
    estDp <- estD / Dmin
  }
  
  n <-  sum(n3x3)
  
  corr <- estD / sqrt( pA * pB * pa * pb )
  
  dchi <- (2*n*estD^2)/(pA * pa * pB* pb)
  dpval <- 1 - pchisq(dchi,1)
  
  retval <- list(
    "D"=estD,
    "D'"=estDp,
    "r" = corr,
    "R^2" = corr^2,
    "n"=n,
    "X^2"=dchi,
    "P-value"=dpval
  )
  retval
}