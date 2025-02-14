###############################
## Genomic Selection
## k-fold cross-validation
###############################

rm(list=ls()) 
setwd("C:/Users/sgeza/OneDrive/Desktop/Online_Training/GS_GWAS/Datasets/SalmonGIll/GSRR")
library(ASRgenomics)
library(BGLR)

# Reading Phenotypic data
data(pheno.salmon)
datag <- pheno.salmon
head(datag)

# Reading full marker data
markerFULL <- read.table('../Mclean.txt', header=TRUE)
markerFULL[1:5,1:5]
dim(markerFULL)  # 1481 x 11146
# NOTE: order does not match 
markerFULL <- markerFULL[order(rownames(markerFULL)),]
markerFULL[1:5,1:5]
sum(rownames(markerFULL) == datag$ID)

##############################

# Fitting a Model - Bayes B

#nIter <- 15000; burnIn <- 5000  # Minimum recommended
nIter <- 6000; burnIn <-2000    # For Testing
X <- as.matrix(markerFULL)
y <- datag$mean_gill_score

ETA <- list(list(X=X, model="BayesB"))

#############################
# k-fold cross-validation

k <- 5
n <- length(datag$ID)
sel <- rep(1:k, length.out = n) 
group <- sample(sel, n)
table(group)  

ypred_cv <- matrix(data=NA, nrow=n, ncol=1)
h2_cv <- matrix(data=NA, nrow=k, ncol=1)

for (g in 1:k)
{ 

  ycv<-y
  
  for (j in 1:n) {
    if(group[j] == g) { ycv[j] <- NA } 
  }
  
  gsBayesBcv <- BGLR(y=ycv, ETA=ETA, nIter=nIter, burnIn=burnIn)
  predGScv <- gsBayesBcv$yHat

  Vary <- var(y, na.rm=TRUE)
  VarE <- gsBayesBcv$varE
  h2_cv[g] <- 1 - VarE/Vary
  
  for (j in 1:n) {
    if(group[j] == g) { ypred_cv[j] <- predGScv[j] } 
  }
}

(corr_cv <- cor(y, ypred_cv, method='pearson', use="complete.obs"))
plot(y, ypred_cv)
h2_cv
mean(h2_cv)

# Predictive Ability  corr(yadj,ghat)
(PA <- corr_cv)

# Predictive Accuracy corr(greal,ghat)
(ACC <- PA/sqrt(mean(h2_cv)))  # Approximation
